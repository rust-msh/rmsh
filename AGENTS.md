# AGENTS.md

## Python tooling

使用 uv 管理 Python。所有 Python 依赖、虚拟环境、构建和发布都通过 `uv` 完成。

- 虚拟环境：项目根目录下使用 `uv venv` 创建
- 依赖管理：`pyproject.toml` (如有) 或 `uv add` 添加依赖
- 构建 Python 包：`uv build`
- 发布：`uv publish`
