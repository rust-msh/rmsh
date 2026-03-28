# emstudio
1. 它是一个基于egui+eframe的桌面应用，同时能编译成wasm运行在浏览器中，有trunk启动
2. 它有standalone模式和cloud模式，默认standalone模式运行
3. 它用于电磁仿真的建模，求解，和结果展示。ui风格是四周是dockable面板(多面板可以tab页切换)，中央区域是tab页可切换，可横向纵向分裂的view panel
4. 窗口顶部菜单栏下面的工具栏采用ribbon bar，可以参考https://github.com/gnibuoz/QRibbon
5. 需要3d模型展示，使用wgpu来渲染
6. 采用cargo workspace形式，自定义组件放到可components crate中
