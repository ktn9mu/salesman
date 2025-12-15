Starter code and data for traveling salesman problem


Files in this directory:

* datareader.cpp : example code to read in the data files (use Makefile)
* datareader.py  : example code to read in the data files
* cities23.dat : list of coordinates for 23 cities in North America
* cities150.dat : 150 cities in North America
* cities1k.dat : 1207 cities in North America
* cities2k.dat : 2063 cities around the world
* routeplot.py : code to plot the globe and salesman's path<br>
usage:<br>
python routeplot.py cities.dat [cities2.dat] -r [="NA"],"World"'<br>
NA = North America, World = Mercator projection of the whole earth
* earth.C : (just for fun) plotting the globe in ROOT


Work in progress:

cities150 317298.645km 49501.533km 28.314s
cities1k 732177.737km 111640.589km 255.087s
cities2k (not yet)

How to run: 
cities150:    ./sales cities150.dat route150.dat an150.csv
cities1k:    ./sales cities1k.dat route1k.dat an1k.csv