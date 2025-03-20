#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File  : png_tif.py
# Author: zdy
# Date  : 2025/3/20
import json
import os

import numpy as np
import rasterio
from PIL import Image

from rasterio.transform import from_bounds

path = r"C:\Users\zdy\Desktop\map2"
with open(os.path.join(path, r'output.json'), "r", encoding="utf-8") as f:
    dict_data = json.load(f)  # 确保它是 Python 字典/列表
    for key, value in dict_data.items():
        z = 19
        x, y = key.split("_")
        # x, y = int(x), int(y)
        img_path = path + "\\{z}\\{x}\\".format(z=z, x=x)
        input_file = img_path + "\\{y}.png".format(y=y)
        output_file = path + '\\tif' + "\\{x}_{y}.tif".format(x=x, y=y)
        with Image.open(input_file) as img:
            img_array = np.array(img)

        lt_lon, rb_lat, rb_lon, lt_lat = value[0][0], value[1][1], value[1][0], value[0][1]
        # lt_lon, rb_lat, rb_lon, lt_lat = value[1][0], value[1][1], value[0][0], value[0][1]
        transform = from_bounds(lt_lon, rb_lat, rb_lon, lt_lat, img_array.shape[1], img_array.shape[0])

        with rasterio.open(
                output_file, 'w',
                driver='GTiff',
                height=img_array.shape[0],
                width=img_array.shape[1],
                count=3,  # RGB 三个波段
                dtype=img_array.dtype,
                crs="EPSG:4326",  # WGS 84 坐标
                transform=transform
        ) as dst:
            dst.write(img_array[:, :, 0], 1)  # 红色波段
            dst.write(img_array[:, :, 1], 2)  # 绿色波段
            dst.write(img_array[:, :, 2], 3)  # 蓝色波段
