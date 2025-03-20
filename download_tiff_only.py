import json
import math
import os
import requests
import cv2
import numpy as np
import geopandas as gpd
import mercantile  # 处理瓦片坐标
from PIL import Image
from shapely.geometry import LineString, box, Polygon
from tqdm import tqdm

from download_tile import myThread


def swap(a, b):
    a, b = b, a
    return a, b


class Point:
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat


class Box:
    def __init__(self, point_lt, point_rb):
        self.point_lt = point_lt
        self.point_rb = point_rb


def build_url(x, y, z):
    return "https://gac-geo.googlecnapps.club/maps/vt?lyrs=s&x={x}&y={y}&z={z}".format(x=x, y=y, z=z)


def download(x, y, z, path):
    proxies = {
        "http": "http://127.0.0.1:7890",
        "https": "http://127.0.0.1:7890"
    }
    url = build_url(x, y, z)
    print(url)
    path = path + "\\{z}\\{x}\\".format(z=z, x=x)
    if not os.path.exists(path):
        os.makedirs(path)
    filepath = path + "\\{y}.png".format(y=y)
    if os.path.exists(filepath) and os.path.getsize(filepath) > 400:
        # print(os.path.getsize(filepath))
        print("skip")
        pass
    else:
        for x in range(0,3):
            response = requests.get(url, proxies=proxies)
            if response.status_code == 200:
                with open(filepath, "wb") as f:
                    f.write(response.content)
                break;
            else:
                print(response.text)
                print("network error!")


def xyz2lonlat(x, y, z):
    n = math.pow(2, z)
    lon = x / n * 360.0 - 180.0
    lat = math.atan(math.sinh(math.pi * (1 - 2 * y / n)))
    lat = lat * 180.0 / math.pi
    return lon, lat


def lonlat2xyz(lon, lat, zoom):
    n = math.pow(2, zoom)
    x = ((lon + 180) / 360) * n
    y = (1 - (math.log(math.tan(math.radians(lat)) + (1 / math.cos(math.radians(lat)))) / math.pi)) / 2 * n
    return int(x), int(y)


def cal_tiff_box(x1, y1, x2, y2, z):
    LT = xyz2lonlat(x1, y1, z)
    RB = xyz2lonlat(x2 + 1, y2 + 1, z)
    return Point(LT[0], LT[1]), Point(RB[0], RB[1])


def cal_tiff_box_from_filtered_tiles(tile_set, z):
    """
    计算基于长城矢量筛选出的瓦片集合的地理边界
    :param tile_set: 瓦片集合 [(x, y, z), ...]
    :param z: 瓦片缩放级别
    :return: (左上角经纬度, 右下角经纬度)
    """
    x_vals = [x for x, y, _ in tile_set]
    y_vals = [y for x, y, _ in tile_set]

    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)

    # 计算左上角、右下角的经纬度
    LT = xyz2lonlat(x_min, y_min, z)  # 左上角 (Longitude, Latitude)
    RB = xyz2lonlat(x_max + 1, y_max + 1, z)  # 右下角 (Longitude, Latitude) +1 确保完整覆盖

    return Point(LT[0], LT[1]), Point(RB[0], RB[1])

def tile_to_lonlat(x, y, z):
    """
    将瓦片编号 (x, y, z) 转换为经纬度 (lon, lat)
    :param x: 瓦片 X 坐标
    :param y: 瓦片 Y 坐标
    :param z: 瓦片 Zoom 级别
    :return: (lon, lat)
    """
    # 瓦片总数
    n = 2.0 ** z

    # 计算经度 (Longitude)
    lon = x / n * 360.0 - 180.0

    # 计算纬度 (Latitude)
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * y / n)))
    lat = math.degrees(lat_rad)

    return lon, lat

def merge(x1, y1, x2, y2, z, path):
    row_list = list()
    for i in range(x1, x2+1):
        col_list = list()
        for j in range(y1, y2+1):
            img_path = path + "\\{z}\\{i}\\{j}.png".format(i=i, j=j, z=z)
            if os.path.exists(img_path):
                col_list.append(cv2.imread(img_path))
            else:
                temp = np.zeros((256,256,3), dtype=np.uint8)
                col_list.append(temp)
        k = np.vstack(col_list)
        row_list.append(k)
    result = np.hstack(row_list)
    cv2.imwrite(path + "//merge.png", result)

def core(path, shp_parh, z):
    # 读取长城矢量数据（Shapefile / GeoJSON）
    gdf = gpd.read_file(shp_parh)  # 确保是 EPSG:4326 (WGS 84)

    # 获取矢量的边界框 (Bounding Box)
    minx, miny, maxx, maxy = gdf.total_bounds  # 获取长城的整体范围
    count = (maxx-minx+1) * (maxy-miny+1)
    # 计算该范围内的所有可能瓦片
    tile_set = set()
    for tile in mercantile.tiles(minx, miny, maxx, maxy, z):
        tile_set.add((tile.x, tile.y, z))  # 存入瓦片编号 (x, y, z)

    # 生成瓦片边界的 Polygon
    def tile_polygon(x, y, z):
        bounds = mercantile.bounds(x, y, z)  # 获取瓦片边界
        return Polygon([(bounds.west, bounds.south), (bounds.west, bounds.north),
                        (bounds.east, bounds.north), (bounds.east, bounds.south)])

    # 只保留与长城相交的瓦片
    data_dict = {}
    filtered_tiles = set()
    for x, y, z in tqdm(tile_set):
        tile_geom = tile_polygon(x, y, z)  # 转换瓦片为矢量 Polygon
        if gdf.geometry.intersects(tile_geom).any():  # 只选取相交的瓦片
            filtered_tiles.add((x, y, z, path))
            LT = xyz2lonlat(x, y, z)  # 左上角
            RB = xyz2lonlat(x + 1, y + 1, z)  # 右下角（确保完整覆盖）
            data_dict[str(x) + '_' + str(y)] = [(LT[0], LT[1]), (RB[0], RB[1])]
    print(len(filtered_tiles))

    # 下载符合条件的瓦片
    # for x, y, z, _ in filtered_tiles:
    #     LT = xyz2lonlat(x, y, z)  # 左上角
    #     RB = xyz2lonlat(x + 1, y + 1, z)  # 右下角（确保完整覆盖）
    #     data_dict[str(x)+'_'+str(y)] = [(LT[0], LT[1]), (RB[0], RB[1])]
        # download(x, y, z, path)
        # print("{m}/{n}".format(m=count, n=all))
        # pass
    json_output = json.dumps(data_dict, indent=4, ensure_ascii=False)
    with open(os.path.join(path, 'output.json'), "w", encoding="utf-8") as f:
        f.write(json_output)

    filteredArray = list(filtered_tiles)
    urlArraySplit = np.array_split(np.array(filteredArray), 16)

    threads = []

    for item in urlArraySplit:
        thread = myThread(item)
        thread.start()
        threads.append(thread)

    for thread in threads:
        thread.join()

    # 获取所有 x, y 的范围
    # x_vals = [x for x, y, _ in tile_set]
    # y_vals = [y for x, y, _ in tile_set]
    #
    # x_min, x_max = min(x_vals), max(x_vals)
    # y_min, y_max = min(y_vals), max(y_vals)
    # merge(x_min, y_min, x_max, y_max, z, path)
    #
    # lt, rb = cal_tiff_box_from_filtered_tiles(filtered_tiles, z)
    # cmd = "gdal_translate.exe -of GTiff -a_srs EPSG:4326 -a_ullr {p1_lon} " \
    #       "{p1_lat} {p2_lon} {p2_lat}" \
    #       " {input} {output}".format(p1_lon=lt.lon, p1_lat=lt.lat, p2_lon=rb.lon, p2_lat=rb.lat,
    #                                  input='/'.join(path.split('\\'))+"/merge.png", output='/'.join(path.split('\\'))+"/output.tiff")
    #
    # print('配置环境变量 然后运行 即可生成 tiff \n' + cmd)


if __name__ == '__main__':

    path = r"H:\googlemaps\map_河北唐山_秦皇岛卢龙_北京平谷"
    shp_parh = r"E:\PycharmProjects\spider\长城重点区段矢量\河北唐山_秦皇岛卢龙_北京平谷.shp"

    # 需要生成json文件！！！！！！！！！！！！！！！！！！！！
    core(path, shp_parh, z = 19) #调整下载级别


