# google_map_satellite_download

#0.1 
  上家公司有内部使用的地图下载工具，支持下载各种影像，离职后，一次工作中使用无偏移的谷歌影像
，发现市场上的软件如bigemap之类收费有水印，有些不便，于是自己动手写了一个。简单实现了下面三个功能。

- 支持下载谷歌地图卫星影像散列瓦片
- 支持下载谷歌地图卫星合并为一张png 
- 支持合并后图片转换为geotiff 

#0.2 使用
1. 需要能访问google服务器，在代码中配置你的proxy
2. 如果需要生成geotiff文件需要安装GDAL，[已经编译好的GDAL windows版本](https://www.gisinternals.com/release.php)
3. 如果不需要合并或者生成geotiff可以注释掉相代码

本职是后台开发，亦不是GIS专业，所以对OpenCV，GDAL，Python也没有深入了解，见谅。



