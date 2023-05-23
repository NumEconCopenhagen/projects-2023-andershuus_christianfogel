# Data analysis project

Our project is titled **The correlation between family income and membersship of the Danish church.** and is about examining whether there can be found a correlation between the being a member of the Danish church and average family disposable income.

The **results** of the project can be seen from running [dataproject.ipynb](dataproject.ipynb).

We download the **following datasets** from the DST API:

1: INDKP132 - this has information on disposable family income
2: KM6 - this has information on membership status in the Danish church

We use a package, found on github, that can produce a map of Denmark: https://github.com/Neogeografen/dagi . This is the downloaded file "kommuner.geojson" in the folder. 

**Dependencies:** Apart from a standard Anaconda Python 3 installation, the project requires the following installations:

``%pip install geopandas``
``%pip install mapclassify``'
``%pip install folium``
``%pip install branca``