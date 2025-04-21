# Methodology-for-the-Identification-of-Ground-Water-Flow-Types-in-Caverns-using-LiDAR-and-MATLAB-
Uses a MATLAB script to identify stalactites then analyzes morphological properties of these stalactites to classify them as 1 of 3 flow types based on aspect ration and cross sectional area. Developed by DR. Kashif Mahmude and Tested by Matthew Kaspar and Kathryn Brown. (Repository of the Script and the base Lidar data from Natural Bridge Caverns)

This method takes a LiDAr bin file and inmterpolates it into an .SGEMS file which can then be input into the script to identify stalactites and flow paths.

To identify stalactites, the method creates an averaged surface to represent the overall ceiling topography, using the moving window average of the ceiling surface. We tested different window size to find the optimum values for our site.  Once we have the smoothed ceiling surface, the next step is to subtract the actual ceiling point cloud from the smoothed ceiling surface to generate map of anomalies. These anomalies are then defined the location of stalactites within the ceiling section. We then plotted histogram of the anomalies, and used a threshold to identify the individual stalactites. Anything below this threshold were natural rock surface variability and were not identified as stalactites. Lastly, we performed a connected component analysis (Mahmud et al. 2015) with the flow percentile threshold to connect the nearby stalactites in a group. We then use the morphological properties of those groups of stalactites to determine the different water flow types following Mahmud et al. (2015). These properties consist of stalactite’s cross-sectional area, and its aspect ratio which is the ratio between its major and minor axes.

the flow classifications are as follows
Matrix flow, flower moving water thast is marked by solitary smaller stalactites such as soda straws
Fracture Flow, which is faster moving water and forms larger linear stalactites like cave ribbon
Combination flow, a mixture of matrix and fracture flow which forms more circular stalactites
