//! Lightweight recursive-descent expression parser and evaluator.
//!
//! Supports:
//! - Arithmetic: `+`, `-`, `*`, `/`, unary `-`, parentheses
//! - Functions: `sqrt`, `sin`, `cos`, `tan`, `exp`, `log`, `abs`, `min`, `max`
//! - Variable references: `$freq` (project), `patch_l` (design)
//! - Numeric literals with optional unit suffix: `2.4GHz`, `28.5mm`, `3.14`

use std::collections::HashMap;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum ExprError {
    #[error("unexpected character '{0}' at position {1}")]
    UnexpectedChar(char, usize),
    #[error("unexpected end of expression")]
    UnexpectedEnd,
    #[error("undefined variable: {0}")]
    UndefinedVariable(String),
    #[error("unknown function: {0}")]
    UnknownFunction(String),
    #[error("circular dependency involving: {0}")]
    CircularDependency(String),
    #[error("parse error: {0}")]
    ParseError(String),
}

// ---------------------------------------------------------------------------
// Unit multiplier parsing
// ---------------------------------------------------------------------------

/// Parse a unit suffix and return the SI multiplier.
pub fn unit_multiplier(suffix: &str) -> f64 {
    match suffix {
        // Frequency
        "Hz" => 1.0,
        "kHz" | "KHz" => 1e3,
        "MHz" => 1e6,
        "GHz" => 1e9,
        "THz" => 1e12,
        // Length
        "m" => 1.0,
        "cm" => 1e-2,
        "mm" => 1e-3,
        "um" | "µm" => 1e-6,
        "nm" => 1e-9,
        "mil" => 25.4e-6,
        "in" => 25.4e-3,
        // Angle
        "deg" => std::f64::consts::PI / 180.0,
        "rad" => 1.0,
        // Electrical
        "ohm" | "Ohm" => 1.0,
        "kohm" | "kOhm" => 1e3,
        "Mohm" | "MOhm" => 1e6,
        "S" => 1.0,
        "mS" => 1e-3,
        // Temperature (passthrough — just store the number)
        "cel" | "K" => 1.0,
        // Time
        "s" => 1.0,
        "ms" => 1e-3,
        "us" | "µs" => 1e-6,
        "ns" => 1e-9,
        "ps" => 1e-12,
        // Power / miscellaneous
        "W" => 1.0,
        "mW" => 1e-3,
        "dB" | "dBm" | "dBi" => 1.0, // passthrough for dB values
        // Capacitance
        "F" => 1.0,
        "pF" => 1e-12,
        "nF" => 1e-9,
        "uF" | "µF" => 1e-6,
        // Inductance
        "H" => 1.0,
        "nH" => 1e-9,
        "uH" | "µH" => 1e-6,
        "mH" => 1e-3,
        _ => 1.0,
    }
}

/// Parse a value string like `"2.4GHz"` or `"28.5mm"` into an SI f64.
pub fn parse_value_with_unit(s: &str) -> Result<f64, ExprError> {
    let s = s.trim();
    if s.is_empty() {
        return Err(ExprError::ParseError("empty value".into()));
    }

    // Find where the numeric part ends
    let num_end = s
        .find(|c: char| c.is_alphabetic() || c == 'µ')
        .unwrap_or(s.len());

    let num_str = &s[..num_end];
    let unit_str = &s[num_end..];

    let number: f64 = num_str
        .parse()
        .map_err(|_| ExprError::ParseError(format!("invalid number: {num_str}")))?;

    if unit_str.is_empty() {
        Ok(number)
    } else {
        Ok(number * unit_multiplier(unit_str))
    }
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(f64),
    Ident(String),    // variable or function name
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
    Comma,
    Eof,
}

struct Tokenizer {
    chars: Vec<char>,
    pos: usize,
}

impl Tokenizer {
    fn new(input: &str) -> Self {
        Self {
            chars: input.chars().collect(),
            pos: 0,
        }
    }

    fn peek_char(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<char> {
        let c = self.chars.get(self.pos).copied();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn skip_whitespace(&mut self) {
        while self.peek_char().is_some_and(|c| c.is_whitespace()) {
            self.advance();
        }
    }

    fn next_token(&mut self) -> Result<Token, ExprError> {
        self.skip_whitespace();
        match self.peek_char() {
            None => Ok(Token::Eof),
            Some(c) => match c {
                '+' => { self.advance(); Ok(Token::Plus) }
                '-' => { self.advance(); Ok(Token::Minus) }
                '*' => { self.advance(); Ok(Token::Star) }
                '/' => { self.advance(); Ok(Token::Slash) }
                '(' => { self.advance(); Ok(Token::LParen) }
                ')' => { self.advance(); Ok(Token::RParen) }
                ',' => { self.advance(); Ok(Token::Comma) }
                '$' | '_' if !c.is_ascii_digit() => {
                    // Variable or identifier starting with $ or _
                    self.read_ident()
                }
                _ if c.is_ascii_digit() || c == '.' => self.read_number(),
                _ if c.is_alphabetic() || c == '_' => self.read_ident(),
                _ => Err(ExprError::UnexpectedChar(c, self.pos)),
            },
        }
    }

    fn read_number(&mut self) -> Result<Token, ExprError> {
        let start = self.pos;
        while self.peek_char().is_some_and(|c| c.is_ascii_digit() || c == '.' || c == 'e' || c == 'E') {
            self.advance();
            // Handle exponent sign: e.g. 1e-3
            if self.chars.get(self.pos.wrapping_sub(1)).is_some_and(|&c| c == 'e' || c == 'E') {
                if self.peek_char().is_some_and(|c| c == '+' || c == '-') {
                    self.advance();
                }
            }
        }

        // Check for unit suffix (letters immediately after the number)
        let num_end = self.pos;
        let unit_start = self.pos;
        while self.peek_char().is_some_and(|c| c.is_alphabetic() || c == 'µ') {
            self.advance();
        }

        let num_str: String = self.chars[start..num_end].iter().collect();
        let number: f64 = num_str
            .parse()
            .map_err(|_| ExprError::ParseError(format!("invalid number: {num_str}")))?;

        if self.pos > unit_start {
            let unit: String = self.chars[unit_start..self.pos].iter().collect();
            Ok(Token::Number(number * unit_multiplier(&unit)))
        } else {
            Ok(Token::Number(number))
        }
    }

    fn read_ident(&mut self) -> Result<Token, ExprError> {
        let start = self.pos;
        // Allow $ as first character
        if self.peek_char() == Some('$') {
            self.advance();
        }
        while self.peek_char().is_some_and(|c| c.is_alphanumeric() || c == '_') {
            self.advance();
        }
        let name: String = self.chars[start..self.pos].iter().collect();
        Ok(Token::Ident(name))
    }
}

// ---------------------------------------------------------------------------
// Parser (recursive descent)
// ---------------------------------------------------------------------------

/// Evaluate a mathematical expression with variable substitution.
///
/// # Grammar
/// ```text
/// expr     = term (('+' | '-') term)*
/// term     = unary (('*' | '/') unary)*
/// unary    = '-' unary | primary
/// primary  = NUMBER | IDENT '(' args ')' | IDENT | '(' expr ')'
/// args     = expr (',' expr)*
/// ```
pub fn evaluate(input: &str, vars: &HashMap<String, f64>) -> Result<f64, ExprError> {
    let mut tokenizer = Tokenizer::new(input);
    let mut tokens = Vec::new();
    loop {
        let tok = tokenizer.next_token()?;
        let is_eof = tok == Token::Eof;
        tokens.push(tok);
        if is_eof {
            break;
        }
    }

    let mut parser = Parser { tokens, pos: 0 };
    let result = parser.parse_expr(vars)?;

    if parser.peek() != &Token::Eof {
        return Err(ExprError::ParseError("trailing tokens".into()));
    }

    Ok(result)
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn peek(&self) -> &Token {
        self.tokens.get(self.pos).unwrap_or(&Token::Eof)
    }

    fn advance(&mut self) -> Token {
        let tok = self.tokens.get(self.pos).cloned().unwrap_or(Token::Eof);
        self.pos += 1;
        tok
    }

    fn expect(&mut self, expected: &Token) -> Result<(), ExprError> {
        let tok = self.advance();
        if std::mem::discriminant(&tok) == std::mem::discriminant(expected) {
            Ok(())
        } else {
            Err(ExprError::ParseError(format!(
                "expected {expected:?}, got {tok:?}"
            )))
        }
    }

    fn parse_expr(&mut self, vars: &HashMap<String, f64>) -> Result<f64, ExprError> {
        let mut left = self.parse_term(vars)?;
        loop {
            match self.peek() {
                Token::Plus => {
                    self.advance();
                    left += self.parse_term(vars)?;
                }
                Token::Minus => {
                    self.advance();
                    left -= self.parse_term(vars)?;
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_term(&mut self, vars: &HashMap<String, f64>) -> Result<f64, ExprError> {
        let mut left = self.parse_unary(vars)?;
        loop {
            match self.peek() {
                Token::Star => {
                    self.advance();
                    left *= self.parse_unary(vars)?;
                }
                Token::Slash => {
                    self.advance();
                    left /= self.parse_unary(vars)?;
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_unary(&mut self, vars: &HashMap<String, f64>) -> Result<f64, ExprError> {
        if *self.peek() == Token::Minus {
            self.advance();
            Ok(-self.parse_unary(vars)?)
        } else {
            self.parse_primary(vars)
        }
    }

    fn parse_primary(&mut self, vars: &HashMap<String, f64>) -> Result<f64, ExprError> {
        match self.advance() {
            Token::Number(n) => Ok(n),
            Token::Ident(name) => {
                // Check if it's a function call
                if *self.peek() == Token::LParen {
                    self.advance(); // consume '('
                    let mut args = vec![self.parse_expr(vars)?];
                    while *self.peek() == Token::Comma {
                        self.advance();
                        args.push(self.parse_expr(vars)?);
                    }
                    self.expect(&Token::RParen)?;
                    call_function(&name, &args)
                } else {
                    // Variable lookup
                    vars.get(&name)
                        .copied()
                        .ok_or(ExprError::UndefinedVariable(name))
                }
            }
            Token::LParen => {
                let val = self.parse_expr(vars)?;
                self.expect(&Token::RParen)?;
                Ok(val)
            }
            other => Err(ExprError::ParseError(format!("unexpected token: {other:?}"))),
        }
    }
}

fn call_function(name: &str, args: &[f64]) -> Result<f64, ExprError> {
    match name {
        "sqrt" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].sqrt())
        }
        "sin" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].sin())
        }
        "cos" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].cos())
        }
        "tan" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].tan())
        }
        "exp" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].exp())
        }
        "log" | "ln" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].ln())
        }
        "log10" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].log10())
        }
        "abs" => {
            ensure_arg_count(name, args, 1)?;
            Ok(args[0].abs())
        }
        "min" => {
            if args.len() < 2 {
                return Err(ExprError::ParseError(format!(
                    "min requires at least 2 arguments, got {}",
                    args.len()
                )));
            }
            Ok(args.iter().copied().fold(f64::INFINITY, f64::min))
        }
        "max" => {
            if args.len() < 2 {
                return Err(ExprError::ParseError(format!(
                    "max requires at least 2 arguments, got {}",
                    args.len()
                )));
            }
            Ok(args.iter().copied().fold(f64::NEG_INFINITY, f64::max))
        }
        "pow" => {
            ensure_arg_count(name, args, 2)?;
            Ok(args[0].powf(args[1]))
        }
        _ => Err(ExprError::UnknownFunction(name.to_string())),
    }
}

fn ensure_arg_count(name: &str, args: &[f64], expected: usize) -> Result<(), ExprError> {
    if args.len() != expected {
        Err(ExprError::ParseError(format!(
            "{name} expects {expected} argument(s), got {}",
            args.len()
        )))
    } else {
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Variable resolution with topological sorting
// ---------------------------------------------------------------------------

use crate::variable::Variable;

/// Resolve all variables in dependency order, returning a map of name → SI value.
///
/// Project variables are prefixed with `$`, design variables are not.
pub fn resolve_variables(
    project_vars: &HashMap<String, Variable>,
    design_vars: &HashMap<String, Variable>,
) -> Result<HashMap<String, f64>, ExprError> {
    let mut resolved = HashMap::new();

    // Build a combined list
    let mut all_vars: Vec<(String, &Variable)> = Vec::new();
    for (name, var) in project_vars {
        all_vars.push((name.clone(), var));
    }
    for (name, var) in design_vars {
        all_vars.push((name.clone(), var));
    }

    // Simple iterative resolution: keep trying until no more progress
    let mut remaining = all_vars;
    let mut made_progress = true;

    while made_progress && !remaining.is_empty() {
        made_progress = false;
        let mut still_remaining = Vec::new();

        for (name, var) in remaining {
            if let Some(value_str) = &var.value {
                match parse_value_with_unit(value_str) {
                    Ok(val) => {
                        resolved.insert(name, val);
                        made_progress = true;
                        continue;
                    }
                    Err(_) => {
                        // Try as expression
                        match evaluate(value_str, &resolved) {
                            Ok(val) => {
                                resolved.insert(name, val);
                                made_progress = true;
                                continue;
                            }
                            Err(_) => {
                                still_remaining.push((name, var));
                            }
                        }
                    }
                }
            } else if let Some(expr) = &var.expression {
                match evaluate(expr, &resolved) {
                    Ok(val) => {
                        resolved.insert(name, val);
                        made_progress = true;
                        continue;
                    }
                    Err(_) => {
                        still_remaining.push((name, var));
                    }
                }
            } else {
                // No value or expression — skip (defaults to 0?)
                resolved.insert(name, 0.0);
                made_progress = true;
            }
        }

        remaining = still_remaining;
    }

    if !remaining.is_empty() {
        let names: Vec<_> = remaining.iter().map(|(n, _)| n.as_str()).collect();
        return Err(ExprError::CircularDependency(names.join(", ")));
    }

    Ok(resolved)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_number() {
        let vars = HashMap::new();
        assert!((evaluate("42", &vars).unwrap() - 42.0).abs() < 1e-10);
        assert!((evaluate("3.14", &vars).unwrap() - 3.14).abs() < 1e-10);
    }

    #[test]
    fn parse_number_with_unit() {
        let vars = HashMap::new();
        let val = evaluate("2.4GHz", &vars).unwrap();
        assert!((val - 2.4e9).abs() < 1.0);
        let val = evaluate("28.5mm", &vars).unwrap();
        assert!((val - 28.5e-3).abs() < 1e-10);
    }

    #[test]
    fn arithmetic() {
        let vars = HashMap::new();
        assert!((evaluate("2 + 3", &vars).unwrap() - 5.0).abs() < 1e-10);
        assert!((evaluate("10 - 4 * 2", &vars).unwrap() - 2.0).abs() < 1e-10);
        assert!((evaluate("(10 - 4) * 2", &vars).unwrap() - 12.0).abs() < 1e-10);
        assert!((evaluate("6 / 3 + 1", &vars).unwrap() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn unary_minus() {
        let vars = HashMap::new();
        assert!((evaluate("-5", &vars).unwrap() + 5.0).abs() < 1e-10);
        assert!((evaluate("-(3 + 2)", &vars).unwrap() + 5.0).abs() < 1e-10);
    }

    #[test]
    fn variables() {
        let mut vars = HashMap::new();
        vars.insert("$freq".to_string(), 2.4e9);
        vars.insert("patch_l".to_string(), 28.5e-3);

        let val = evaluate("$freq / 1e9", &vars).unwrap();
        assert!((val - 2.4).abs() < 1e-10);

        let val = evaluate("patch_l * 2", &vars).unwrap();
        assert!((val - 57.0e-3).abs() < 1e-10);
    }

    #[test]
    fn functions() {
        let vars = HashMap::new();
        let val = evaluate("sqrt(4)", &vars).unwrap();
        assert!((val - 2.0).abs() < 1e-10);

        let val = evaluate("abs(-3.5)", &vars).unwrap();
        assert!((val - 3.5).abs() < 1e-10);

        let val = evaluate("min(3, 1, 2)", &vars).unwrap();
        assert!((val - 1.0).abs() < 1e-10);

        let val = evaluate("max(3, 1, 2)", &vars).unwrap();
        assert!((val - 3.0).abs() < 1e-10);

        let val = evaluate("pow(2, 10)", &vars).unwrap();
        assert!((val - 1024.0).abs() < 1e-10);
    }

    #[test]
    fn parse_value_with_unit_cases() {
        assert!((parse_value_with_unit("2.4GHz").unwrap() - 2.4e9).abs() < 1.0);
        assert!((parse_value_with_unit("28.5mm").unwrap() - 28.5e-3).abs() < 1e-10);
        assert!((parse_value_with_unit("50ohm").unwrap() - 50.0).abs() < 1e-10);
        assert!((parse_value_with_unit("22cel").unwrap() - 22.0).abs() < 1e-10);
        assert!((parse_value_with_unit("3.14").unwrap() - 3.14).abs() < 1e-10);
    }

    #[test]
    fn resolve_variables_basic() {
        let mut project = HashMap::new();
        project.insert(
            "$freq".to_string(),
            Variable {
                value: Some("2.4GHz".to_string()),
                expression: None,
                description: String::new(),
                unit_type: Some("Frequency".to_string()),
            },
        );

        let mut design = HashMap::new();
        design.insert(
            "half_freq".to_string(),
            Variable {
                value: None,
                expression: Some("$freq / 2".to_string()),
                description: String::new(),
                unit_type: Some("Frequency".to_string()),
            },
        );

        let resolved = resolve_variables(&project, &design).unwrap();
        assert!((resolved["$freq"] - 2.4e9).abs() < 1.0);
        assert!((resolved["half_freq"] - 1.2e9).abs() < 1.0);
    }

    #[test]
    fn circular_dependency_detected() {
        let mut vars = HashMap::new();
        vars.insert(
            "a".to_string(),
            Variable {
                value: None,
                expression: Some("b + 1".to_string()),
                description: String::new(),
                unit_type: None,
            },
        );
        vars.insert(
            "b".to_string(),
            Variable {
                value: None,
                expression: Some("a + 1".to_string()),
                description: String::new(),
                unit_type: None,
            },
        );

        let result = resolve_variables(&HashMap::new(), &vars);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), ExprError::CircularDependency(_)));
    }
}
