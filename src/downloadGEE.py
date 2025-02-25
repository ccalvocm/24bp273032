import os
import ee
import pandas as pd
import datetime
import requests
import zipfile
import xarray as xr


def time_index_from_filenames(filenames):
    '''helper function to create a pandas DatetimeIndex
       Filename example: 20150520_0164.tif'''
    return pd.DatetimeIndex([pd.Timestamp(f[:8]) for f in filenames])


def files2NETCDF(folderName):
    filenames = os.listdir(folderName)
    files = sorted([x for x in filenames if ('.tif' in x) & ('xml' not in x)])
    time = xr.Variable('time', time_index_from_filenames(files))
    chunks = {'x': 5490, 'y': 5490, 'band': 1}
    da = xr.concat([xr.open_dataset(os.path.join(
        '.', folderName, f), chunks=chunks) for f in files], dim=time)
    da.to_netcdf(os.path.join(folderName, 'netcdf.nc'))


def pythonDate2eeDate(date=datetime.datetime):
    eeDate = date.strftime("%Y-%m-%d")
    return eeDate


def outFolder(folderName):
    try:
        os.mkdir(os.path.join('.', folderName))
    except:
        error = 'directorio ya existe'
    return None

def download(collection, band, folderName):
    # login()

    # 1. authenticate
    ee.Authenticate()
    ee.Initialize()

    # 2. define python dates
    fechaIni = datetime.datetime(1979,4, 1)
    fechaFin = datetime.datetime(2022, 3, 31)
    # print(fechaIni, fechaFin)
    # print(pythonDate2eeDate(fechaFin))

    # 3. define the collection
    # collection = "ECMWF/ERA5_LAND/DAILY_AGGR"
    imageCollection = ee.ImageCollection(collection)

    # 4. define the bounding box and filter accordingly
    llx, lly, urx, ury = -71.71782, \
    -32.28247,	-69.809361	, -29.0366
    rectangle = ee.Geometry.Rectangle([llx, lly, urx, ury])
    imageCollectionRegion = imageCollection.filterBounds(rectangle)

    # 5. select the band
    # band = "evaporation_from_open_water_surfaces_excluding_oceans_sum"
    imageCollectionBand = imageCollectionRegion.select(band)
    f = open("log.txt", "w")

    # 6. output folder
    outFolder(folderName)
    rutaDl = os.path.join(folderName)

    # 7. select date
    for date in pd.date_range(fechaIni, fechaFin):
        date = pythonDate2eeDate(date)
        print(f"Attempting download image for {date}")
        try:
            imageCollectionDate = imageCollectionBand.filterDate(date).min().clip(rectangle)

            url = imageCollectionDate.getDownloadURL({
                'scale': 11132,
                'crs': 'EPSG:32719',
                'format': 'GEO_TIFF',
                'region': rectangle})
            print(url)
            # print(f"Attempting download image for {date} OK!")
            response = requests.get(url, stream=True)
            filePath = os.path.join(rutaDl, f"tmax{date}.tif")

            # zip
            with open(os.path.join(rutaDl, band+'_'+date+'.tif'), "wb") as fd:
                for chunk in response.iter_content(chunk_size=1024):
                    fd.write(chunk)
            fd.close()

        except:
            # print(f"Attempting download image for {date} FAILED!")
            f.write(f"{date}\n")

    # 8. extract files
    # files = os.listdir(rutaDl)
    # files = [x for x in files if x.endswith('.zip')]
    # for file in files:
    #     zip_ref = zipfile.ZipFile(os.path.join(
    #         rutaDl, file))  # create zipfile object
    #     zip_ref.extractall(rutaDl)  # extract file to dir
    #     zip_ref.close()  # close file
    #     os.remove(os.path.join(rutaDl, file))  # delete zipped file



def outFolder(folderName):
    if not os.path.exists(folderName):
        os.makedirs(folderName)

def download_annual_max(collection, band, folderName, start_year=1979, end_year=2022):
    """Download annual maximum values for ERA5 Land daily aggregated data"""
    
    ee.Authenticate()
    ee.Initialize()
    # Open image collection and select band
    imageCollection = ee.ImageCollection(collection)
    
    # Define region (adjust coordinates as needed)
    llx, lly, urx, ury = -71.71782, -32.28247, -69.809361, -29.0366
    rectangle = ee.Geometry.Rectangle([llx, lly, urx, ury])
    
    imageCollectionFiltered = imageCollection.select(band).filterBounds(rectangle)
    
    outFolder(folderName)
    
    for year in range(start_year, end_year + 1):
        try:
            print(f"Processing year {year}...")
            start_date = ee.Date.fromYMD(year, 1, 1)
            end_date = ee.Date.fromYMD(year + 1, 1, 1)
            
            # Compute annual maximum image
            yearly_max = imageCollectionFiltered.filterDate(start_date, end_date).max().clip(rectangle)
            
            # Get direct download URL. Adjust 'scale' as needed; here 11132 and EPSG:32719 are used.
            url = yearly_max.getDownloadURL({
                'scale': 11132,
                'crs': 'EPSG:32719',
                'format': 'GEO_TIFF',
                'region': rectangle
            })
            
            # Download the file using requests
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            output_file = os.path.join(folderName, f"{band}_max_{year}.tif")
            with open(output_file, "wb") as fd:
                for chunk in response.iter_content(chunk_size=1024*1024):
                    fd.write(chunk)
            
            print(f"Successfully downloaded maximum for year {year}")
        except Exception as e:
            print(f"Error processing year {year}: {e}")
            error_log = os.path.join(folderName, "error_log.txt")
            with open(error_log, "a") as f:
                f.write(f"{year}: {str(e)}\n")

if __name__ == '__main__':
    collection = "ECMWF/ERA5_LAND/DAILY_AGGR"
    band = "runoff_sum"  # Change this to your desired band name
    folderName = os.path.join("/", "Users", "farrospide", "Downloads", "runoff_ERA5_UTM")
    download_annual_max(collection=collection, band=band, folderName=folderName, start_year=1979, end_year=2022)