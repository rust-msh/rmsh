# emstudio
EmStudio是一个类似Ansys electronics desktop的电磁仿真工具，功能对标 Q3D+HFSS, 它后端依赖Rem(https://github.com/javagg/rem2.git)来完成电磁仿真。它采用rust语言开发，采用wgpu进行3d渲染，能在本地和浏览器（通过wasm技术）中运行。

## Web Local-First 验证

执行以下脚本可快速检查 Web 部署链路（WASM 编译、Worker、工具链状态）：

```bash
chmod +x scripts/verify-web-localfirst.sh
./scripts/verify-web-localfirst.sh
```

如需完整产物构建：

```bash
chmod +x scripts/build-wasm.sh
./scripts/build-wasm.sh
```

说明：

- 脚本会分别检查 `emstudio-render`、`emstudio-worker`、`emstudio-main` 的 `wasm32-unknown-unknown` 编译。
- 若 `emstudio-main` 出现 `getrandom` 的 wasm 提示（需 `js` feature），脚本会给出告警而不阻塞其它检查。
