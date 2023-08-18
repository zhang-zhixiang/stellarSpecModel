import os
import time
from hashlib import md5
from . import config


def fetch_grid(grid_name):
    if grid_name not in config.grid_names:
        raise ValueError(f'grid_name should be one of {list(config.grid_names.keys())}')
    file_name, url, md5_value = config.grid_names[grid_name]
    download_dir = config.grid_data_dir
    expected_file_name = file_name
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    if not os.path.exists(os.path.join(download_dir, expected_file_name)):
        download_from_jianguoyun(url, download_dir, expected_file_name)
    if not check_md5(os.path.join(download_dir, expected_file_name), md5_value):
        raise ValueError(f'{expected_file_name} is broken, please delete it and try again')
    return os.path.join(download_dir, expected_file_name)


def check_md5(file_path, md5_value):
    with open(file_path, 'rb') as f:
        md5_obj = md5()
        md5_obj.update(f.read())
        md5_file = md5_obj.hexdigest()
    return md5_file == md5_value


def download_from_jianguoyun(url, download_dir, expected_file_name):
    try:
        import selenium
        from selenium import webdriver
        from selenium.webdriver.chrome.options import Options
    except:
        print(('Please install selenium first: pip install selenium\n'
            'If pip report "error: externally-managed-environment"\n'
            'please use: pip install selenium --user --break-system-packages'))
        raise ImportError('selenium is not installed')
    if selenium.__version__ < '4.11.0':
        raise ImportError('selenium version should be >= 4.11.0')
    # 设置Chrome浏览器下载文件的默认目录（在主目录下的指定文件夹）
    # download_dir = os.path.expanduser("~")

    # 配置Chrome选项
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # 无头模式，不显示浏览器窗口
    chrome_options.add_argument("--disable-popup-blocking")  # 禁用弹出窗口阻止
    chrome_options.add_argument("--disable-notifications")  # 禁用通知
    chrome_options.add_argument("--disable-gpu")  # 禁用GPU加速

    # 设置下载文件的相关配置
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "download.directory_upgrade": True,
        "safebrowsing.enabled": False,
    }

    chrome_options.add_experimental_option("prefs", prefs)

    # 创建Chrome浏览器实例
    driver = webdriver.Chrome(options=chrome_options)

    # 打开一个网页，触发文件下载操作
    # url = 'https://www.jianguoyun.com/p/DZmcNoUQ2ZfcCBjW-5cFIAA'

    driver.get(url)

    # 点击下载链接
    try:
        driver.find_element('partial link text', 'Download').click()
    except:
        driver.find_element('partial link text', '下载').click()

    # 等待一段时间，留足够的时间给文件下载
    t0 = time.time()
    while time.time() - t0 < 60:
        time.sleep(5)  # 这里可以根据实际情况调整等待时间
        print('.', end="", flush=True)
        # 检查下载目录中是否有预期的文件
        # expected_file_name = "MARCS_grid.hdf5"  # 你期望的文件名
        downloaded_files = os.listdir(download_dir)
        file_downloaded = expected_file_name in downloaded_files

        # 如果文件已经下载，则退出循环
        if file_downloaded:
            print(f"File {expected_file_name} downloaded successfully.")
            break

    # 关闭浏览器实例
    driver.quit()

    # 如果文件没有下载，则抛出异常
    if not file_downloaded:
        print(f"File {expected_file_name} not found in download directory.")
        raise FileNotFoundError(f"File {expected_file_name} not found in {download_dir} directory.")
