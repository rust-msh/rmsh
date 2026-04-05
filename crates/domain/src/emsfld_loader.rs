// ---------------------------------------------------------------------------
// EMStudio Field Data (.emsfld) Binary Loader
// ---------------------------------------------------------------------------
//
// Binary format for frequency-domain FEM field solutions.
// Supports random access by frequency index via an index table.
//
// File layout:
//   [Header 128 bytes] [Frequency Table] [Index Table] [Field Blocks...]
//
// See docs/em-result-file-formats.md §3.3 for full specification.

use std::io::{self, Read, Seek, SeekFrom};
use std::path::Path;

use thiserror::Error;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug, Error)]
pub enum FieldError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("Invalid magic: expected EMSFLD")]
    InvalidMagic,
    #[error("Invalid byte order marker")]
    InvalidByteOrder,
    #[error("Unsupported version: {0}")]
    UnsupportedVersion(u32),
    #[error("Unsupported data type: {0}")]
    UnsupportedDataType(u32),
    #[error("Frequency index {index} out of range (max {max})")]
    FrequencyOutOfRange { index: usize, max: usize },
    #[error("Invalid field data: {0}")]
    InvalidData(String),
}

// ---------------------------------------------------------------------------
// Header (128 bytes, little-endian)
// ---------------------------------------------------------------------------

/// Magic bytes for .emsfld files.
pub const EMSFLD_MAGIC: &[u8; 8] = b"EMSFLD\0\0";
/// Byte order marker (little-endian check).
pub const BYTE_ORDER_MARKER: u32 = 0x01020304;

#[derive(Debug, Clone, Copy)]
pub struct EmsFldHeader {
    pub magic: [u8; 8],
    pub version: u32,
    pub byte_order: u32,
    /// 0=E-field, 1=H-field, 2=J-field, 3=Combined
    pub field_type: u32,
    /// 0=complex f64, 1=complex f32
    pub data_type: u32,
    pub num_nodes: u64,
    /// 3 for vector field (x,y,z), 1 for scalar
    pub num_components: u32,
    pub num_frequencies: u32,
    /// 0=Hz, 1=kHz, 2=MHz, 3=GHz
    pub frequency_unit: u32,
    pub freq_table_offset: u64,
    pub index_offset: u64,
    pub data_offset: u64,
    pub mesh_file: [u8; 32],
    pub _reserved: [u8; 12],
}

impl EmsFldHeader {
    pub fn read_from<R: Read>(reader: &mut R) -> Result<Self, FieldError> {
        let mut buf = [0u8; 128];
        reader.read_exact(&mut buf)?;

        let magic: [u8; 8] = buf[0..8].try_into().unwrap();
        if &magic != EMSFLD_MAGIC {
            return Err(FieldError::InvalidMagic);
        }

        let version = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        if version != 1 {
            return Err(FieldError::UnsupportedVersion(version));
        }

        let byte_order = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        if byte_order != BYTE_ORDER_MARKER {
            return Err(FieldError::InvalidByteOrder);
        }

        let field_type = u32::from_le_bytes(buf[16..20].try_into().unwrap());
        let data_type = u32::from_le_bytes(buf[20..24].try_into().unwrap());
        let num_nodes = u64::from_le_bytes(buf[24..32].try_into().unwrap());
        let num_components = u32::from_le_bytes(buf[32..36].try_into().unwrap());
        let num_frequencies = u32::from_le_bytes(buf[36..40].try_into().unwrap());
        let frequency_unit = u32::from_le_bytes(buf[40..44].try_into().unwrap());
        let freq_table_offset = u64::from_le_bytes(buf[44..52].try_into().unwrap());
        let index_offset = u64::from_le_bytes(buf[52..60].try_into().unwrap());
        let data_offset = u64::from_le_bytes(buf[60..68].try_into().unwrap());

        let mut mesh_file = [0u8; 32];
        mesh_file.copy_from_slice(&buf[68..100]);

        let mut _reserved = [0u8; 12];
        _reserved.copy_from_slice(&buf[100..112]);

        Ok(EmsFldHeader {
            magic,
            version,
            byte_order,
            field_type,
            data_type,
            num_nodes,
            num_components,
            num_frequencies,
            frequency_unit,
            freq_table_offset,
            index_offset,
            data_offset,
            mesh_file,
            _reserved,
        })
    }

    /// Write header to a buffer.
    pub fn to_bytes(&self) -> [u8; 128] {
        let mut buf = [0u8; 128];
        buf[0..8].copy_from_slice(&self.magic);
        buf[8..12].copy_from_slice(&self.version.to_le_bytes());
        buf[12..16].copy_from_slice(&self.byte_order.to_le_bytes());
        buf[16..20].copy_from_slice(&self.field_type.to_le_bytes());
        buf[20..24].copy_from_slice(&self.data_type.to_le_bytes());
        buf[24..32].copy_from_slice(&self.num_nodes.to_le_bytes());
        buf[32..36].copy_from_slice(&self.num_components.to_le_bytes());
        buf[36..40].copy_from_slice(&self.num_frequencies.to_le_bytes());
        buf[40..44].copy_from_slice(&self.frequency_unit.to_le_bytes());
        buf[44..52].copy_from_slice(&self.freq_table_offset.to_le_bytes());
        buf[52..60].copy_from_slice(&self.index_offset.to_le_bytes());
        buf[60..68].copy_from_slice(&self.data_offset.to_le_bytes());
        buf[68..100].copy_from_slice(&self.mesh_file);
        buf[100..112].copy_from_slice(&self._reserved);
        buf
    }

    /// Get associated mesh file name as string.
    pub fn mesh_filename(&self) -> String {
        let end = self.mesh_file.iter().position(|&b| b == 0).unwrap_or(32);
        String::from_utf8_lossy(&self.mesh_file[..end]).to_string()
    }

    /// Frequency unit multiplier to Hz.
    pub fn freq_multiplier(&self) -> f64 {
        match self.frequency_unit {
            0 => 1.0,      // Hz
            1 => 1e3,      // kHz
            2 => 1e6,      // MHz
            3 => 1e9,      // GHz
            _ => 1.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Field Block Index (16 bytes per frequency)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy)]
pub struct FieldBlockInfo {
    pub offset: u64,
    pub size_bytes: u64,
}

// ---------------------------------------------------------------------------
// Complex field data
// ---------------------------------------------------------------------------

/// A single complex component value.
#[derive(Debug, Clone, Copy)]
pub struct ComplexValue {
    pub real: f64,
    pub imag: f64,
}

impl ComplexValue {
    pub fn magnitude(&self) -> f64 {
        (self.real * self.real + self.imag * self.imag).sqrt()
    }

    pub fn phase_rad(&self) -> f64 {
        self.imag.atan2(self.real)
    }

    pub fn phase_deg(&self) -> f64 {
        self.phase_rad().to_degrees()
    }
}

/// Field data for a single frequency point.
#[derive(Debug, Clone)]
pub struct FieldBlock {
    pub frequency: f64,
    /// Per-node, per-component complex values.
    /// For a 3-component vector field: data[node_idx * 3 + comp_idx]
    pub data: Vec<ComplexValue>,
    pub num_nodes: usize,
    pub num_components: usize,
}

impl FieldBlock {
    /// Get the complex value at a given node and component.
    pub fn value_at(&self, node_idx: usize, component: usize) -> Option<&ComplexValue> {
        if component >= self.num_components || node_idx >= self.num_nodes {
            return None;
        }
        self.data.get(node_idx * self.num_components + component)
    }

    /// Get the vector magnitude at a given node (for 3-component fields).
    pub fn vector_magnitude(&self, node_idx: usize) -> Option<f64> {
        if self.num_components < 3 || node_idx >= self.num_nodes {
            return None;
        }
        let base = node_idx * self.num_components;
        let vx = self.data[base].magnitude();
        let vy = self.data[base + 1].magnitude();
        let vz = self.data[base + 2].magnitude();
        Some((vx * vx + vy * vy + vz * vz).sqrt())
    }

    /// Get field magnitudes for all nodes (scalar: component magnitude, vector: vector magnitude).
    pub fn magnitudes(&self) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.num_nodes);
        for i in 0..self.num_nodes {
            if self.num_components == 1 {
                result.push(self.data[i].magnitude());
            } else {
                result.push(self.vector_magnitude(i).unwrap_or(0.0));
            }
        }
        result
    }
}

// ---------------------------------------------------------------------------
// EmsFldFile — lazy-loading field data file
// ---------------------------------------------------------------------------

/// Handle for a .emsfld file. Reads header and index eagerly, field blocks lazily.
#[derive(Debug, Clone)]
pub struct EmsFldFile {
    pub header: EmsFldHeader,
    pub frequencies: Vec<f64>,
    pub block_index: Vec<FieldBlockInfo>,
    file_path: std::path::PathBuf,
}

impl EmsFldFile {
    /// Open a .emsfld file and read header + frequency table + index.
    pub fn open(path: &Path) -> Result<Self, FieldError> {
        let mut file = std::fs::File::open(path)?;

        let header = EmsFldHeader::read_from(&mut file)?;

        // Read frequency table
        file.seek(SeekFrom::Start(header.freq_table_offset))?;
        let nf = header.num_frequencies as usize;
        let mut frequencies = Vec::with_capacity(nf);
        for _ in 0..nf {
            let mut buf = [0u8; 8];
            file.read_exact(&mut buf)?;
            frequencies.push(f64::from_le_bytes(buf));
        }

        // Read block index
        file.seek(SeekFrom::Start(header.index_offset))?;
        let mut block_index = Vec::with_capacity(nf);
        for _ in 0..nf {
            let mut buf = [0u8; 16];
            file.read_exact(&mut buf)?;
            block_index.push(FieldBlockInfo {
                offset: u64::from_le_bytes(buf[0..8].try_into().unwrap()),
                size_bytes: u64::from_le_bytes(buf[8..16].try_into().unwrap()),
            });
        }

        Ok(EmsFldFile {
            header,
            frequencies,
            block_index,
            file_path: path.to_path_buf(),
        })
    }

    /// Number of frequency points.
    pub fn num_frequencies(&self) -> usize {
        self.frequencies.len()
    }

    /// Get frequencies in Hz.
    pub fn frequencies_hz(&self) -> Vec<f64> {
        let mult = self.header.freq_multiplier();
        self.frequencies.iter().map(|&f| f * mult).collect()
    }

    /// Load field data for a specific frequency index (random access).
    pub fn load_block(&self, freq_idx: usize) -> Result<FieldBlock, FieldError> {
        if freq_idx >= self.block_index.len() {
            return Err(FieldError::FrequencyOutOfRange {
                index: freq_idx,
                max: self.block_index.len().saturating_sub(1),
            });
        }

        let info = &self.block_index[freq_idx];
        let mut file = std::fs::File::open(&self.file_path)?;
        file.seek(SeekFrom::Start(info.offset))?;

        let num_nodes = self.header.num_nodes as usize;
        let num_components = self.header.num_components as usize;
        let total_values = num_nodes * num_components;

        let data = match self.header.data_type {
            0 => {
                // complex f64: 2 × f64 per value = 16 bytes
                let mut buf = vec![0u8; total_values * 16];
                file.read_exact(&mut buf)?;
                let mut values = Vec::with_capacity(total_values);
                for i in 0..total_values {
                    let offset = i * 16;
                    let real = f64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap());
                    let imag =
                        f64::from_le_bytes(buf[offset + 8..offset + 16].try_into().unwrap());
                    values.push(ComplexValue { real, imag });
                }
                values
            }
            1 => {
                // complex f32: 2 × f32 per value = 8 bytes
                let mut buf = vec![0u8; total_values * 8];
                file.read_exact(&mut buf)?;
                let mut values = Vec::with_capacity(total_values);
                for i in 0..total_values {
                    let offset = i * 8;
                    let real =
                        f32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap()) as f64;
                    let imag =
                        f32::from_le_bytes(buf[offset + 4..offset + 8].try_into().unwrap()) as f64;
                    values.push(ComplexValue { real, imag });
                }
                values
            }
            dt => return Err(FieldError::UnsupportedDataType(dt)),
        };

        Ok(FieldBlock {
            frequency: self.frequencies[freq_idx],
            data,
            num_nodes,
            num_components,
        })
    }
}

// ---------------------------------------------------------------------------
// Builder — create .emsfld files (for testing and solver output)
// ---------------------------------------------------------------------------

/// Build a .emsfld file from field data.
pub fn write_emsfld(
    path: &Path,
    field_type: u32,
    num_components: u32,
    frequency_unit: u32,
    mesh_filename: &str,
    blocks: &[(f64, &[ComplexValue])], // (frequency, data)
    use_f32: bool,
) -> Result<(), FieldError> {
    use std::io::Write;

    let num_frequencies = blocks.len() as u32;
    let num_nodes = if blocks.is_empty() {
        0u64
    } else {
        (blocks[0].1.len() / num_components as usize) as u64
    };

    // Compute offsets
    let header_size = 128u64;
    let freq_table_size = num_frequencies as u64 * 8;
    let index_size = num_frequencies as u64 * 16;
    let freq_table_offset = header_size;
    let index_offset = freq_table_offset + freq_table_size;
    let data_offset = index_offset + index_size;

    let bytes_per_value: u64 = if use_f32 { 8 } else { 16 };
    let block_size = num_nodes * num_components as u64 * bytes_per_value;

    let mut mesh_file = [0u8; 32];
    let name_bytes = mesh_filename.as_bytes();
    let copy_len = name_bytes.len().min(31);
    mesh_file[..copy_len].copy_from_slice(&name_bytes[..copy_len]);

    let header = EmsFldHeader {
        magic: *EMSFLD_MAGIC,
        version: 1,
        byte_order: BYTE_ORDER_MARKER,
        field_type,
        data_type: if use_f32 { 1 } else { 0 },
        num_nodes,
        num_components,
        num_frequencies,
        frequency_unit,
        freq_table_offset,
        index_offset,
        data_offset,
        mesh_file,
        _reserved: [0u8; 12],
    };

    let mut file = std::fs::File::create(path)?;

    // Write header
    file.write_all(&header.to_bytes())?;

    // Pad header to 128 bytes (to_bytes already returns 128 bytes, but
    // the header struct only uses 112 bytes of meaningful data; rest is reserved)
    // Actually to_bytes returns exactly 128 bytes, but we only wrote 112 bytes of data
    // Let's pad to 128
    let header_bytes = header.to_bytes();
    // We already have 128 bytes from to_bytes, but we need to also write the remaining 16 bytes
    // that aren't in the struct. Actually the struct uses: 8+4+4+4+4+8+4+4+4+8+8+8+32+12 = 112
    // But to_bytes writes 128. Let's just seek to make sure.
    file.seek(SeekFrom::Start(header_size))?;

    // Write frequency table
    for &(freq, _) in blocks {
        file.write_all(&freq.to_le_bytes())?;
    }

    // Write index table
    let mut current_offset = data_offset;
    for _ in blocks {
        file.write_all(&current_offset.to_le_bytes())?;
        file.write_all(&block_size.to_le_bytes())?;
        current_offset += block_size;
    }

    // Write field blocks
    for &(_, data) in blocks {
        if use_f32 {
            for cv in data {
                file.write_all(&(cv.real as f32).to_le_bytes())?;
                file.write_all(&(cv.imag as f32).to_le_bytes())?;
            }
        } else {
            for cv in data {
                file.write_all(&cv.real.to_le_bytes())?;
                file.write_all(&cv.imag.to_le_bytes())?;
            }
        }
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Field type names
// ---------------------------------------------------------------------------

pub fn field_type_name(field_type: u32) -> &'static str {
    match field_type {
        0 => "E-field",
        1 => "H-field",
        2 => "J-field",
        3 => "Combined",
        _ => "Unknown",
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_roundtrip() {
        let mut mesh_file = [0u8; 32];
        mesh_file[..10].copy_from_slice(b"final.msh\0");

        let header = EmsFldHeader {
            magic: *EMSFLD_MAGIC,
            version: 1,
            byte_order: BYTE_ORDER_MARKER,
            field_type: 0,
            data_type: 0,
            num_nodes: 100,
            num_components: 3,
            num_frequencies: 5,
            frequency_unit: 3, // GHz
            freq_table_offset: 128,
            index_offset: 168,
            data_offset: 248,
            mesh_file,
            _reserved: [0u8; 12],
        };

        let bytes = header.to_bytes();
        assert_eq!(bytes.len(), 128);

        let mut cursor = std::io::Cursor::new(&bytes[..]);
        let parsed = EmsFldHeader::read_from(&mut cursor).unwrap();

        assert_eq!(parsed.version, 1);
        assert_eq!(parsed.num_nodes, 100);
        assert_eq!(parsed.num_components, 3);
        assert_eq!(parsed.num_frequencies, 5);
        assert_eq!(parsed.frequency_unit, 3);
        assert_eq!(parsed.mesh_filename(), "final.msh");
        assert!((parsed.freq_multiplier() - 1e9).abs() < 1.0);
    }

    #[test]
    fn complex_value_operations() {
        let cv = ComplexValue {
            real: 3.0,
            imag: 4.0,
        };
        assert!((cv.magnitude() - 5.0).abs() < 1e-10);
        assert!((cv.phase_deg() - 53.13010235).abs() < 1e-4);
    }

    #[test]
    fn write_and_read_f64() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.emsfld");

        // Create test data: 4 nodes, 3 components (vector field), 2 frequencies
        let num_nodes = 4usize;
        let num_comp = 3usize;

        let data1: Vec<ComplexValue> = (0..num_nodes * num_comp)
            .map(|i| ComplexValue {
                real: i as f64 * 1.0,
                imag: i as f64 * 0.5,
            })
            .collect();

        let data2: Vec<ComplexValue> = (0..num_nodes * num_comp)
            .map(|i| ComplexValue {
                real: i as f64 * 2.0,
                imag: i as f64 * 1.0,
            })
            .collect();

        write_emsfld(
            &path,
            0, // E-field
            3, // vector
            3, // GHz
            "test.msh",
            &[(1.0, &data1), (2.0, &data2)],
            false, // f64
        )
        .unwrap();

        // Read back
        let file = EmsFldFile::open(&path).unwrap();
        assert_eq!(file.num_frequencies(), 2);
        assert_eq!(file.header.num_nodes, 4);
        assert_eq!(file.header.num_components, 3);
        assert!((file.frequencies[0] - 1.0).abs() < 1e-10);
        assert!((file.frequencies[1] - 2.0).abs() < 1e-10);

        let freqs_hz = file.frequencies_hz();
        assert!((freqs_hz[0] - 1e9).abs() < 1.0);

        // Load first block
        let block = file.load_block(0).unwrap();
        assert_eq!(block.num_nodes, 4);
        assert_eq!(block.num_components, 3);
        assert_eq!(block.data.len(), 12);
        assert!((block.data[0].real - 0.0).abs() < 1e-10);
        assert!((block.data[1].real - 1.0).abs() < 1e-10);

        // Test value_at
        let v = block.value_at(1, 0).unwrap(); // node 1, component x
        assert!((v.real - 3.0).abs() < 1e-10);

        // Test vector magnitude
        let mag = block.vector_magnitude(0).unwrap();
        // node 0: components (0+0i, 1+0.5i, 2+1i)
        // magnitudes: 0, sqrt(1.25), sqrt(5)
        let expected = (0.0f64 + 1.25 + 5.0).sqrt();
        assert!((mag - expected).abs() < 1e-6);

        // Load second block
        let block2 = file.load_block(1).unwrap();
        assert!((block2.data[1].real - 2.0).abs() < 1e-10);
    }

    #[test]
    fn write_and_read_f32() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_f32.emsfld");

        let data: Vec<ComplexValue> = (0..6)
            .map(|i| ComplexValue {
                real: i as f64,
                imag: 0.0,
            })
            .collect();

        write_emsfld(&path, 0, 3, 2, "mesh.msh", &[(2.4, &data)], true).unwrap();

        let file = EmsFldFile::open(&path).unwrap();
        assert_eq!(file.header.data_type, 1); // f32
        let block = file.load_block(0).unwrap();
        assert!((block.data[3].real - 3.0).abs() < 1e-5);
    }

    #[test]
    fn frequency_out_of_range() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.emsfld");

        let data: Vec<ComplexValue> = vec![ComplexValue {
            real: 1.0,
            imag: 0.0,
        }];
        write_emsfld(&path, 0, 1, 0, "", &[(1.0, &data)], false).unwrap();

        let file = EmsFldFile::open(&path).unwrap();
        assert!(file.load_block(5).is_err());
    }

    #[test]
    fn field_block_magnitudes() {
        let block = FieldBlock {
            frequency: 1.0,
            data: vec![
                ComplexValue {
                    real: 3.0,
                    imag: 4.0,
                },
                ComplexValue {
                    real: 0.0,
                    imag: 0.0,
                },
                ComplexValue {
                    real: 0.0,
                    imag: 0.0,
                },
            ],
            num_nodes: 1,
            num_components: 3,
        };

        let mags = block.magnitudes();
        assert_eq!(mags.len(), 1);
        // vector magnitude: sqrt(5^2 + 0 + 0) = 5
        assert!((mags[0] - 5.0).abs() < 1e-10);
    }
}
