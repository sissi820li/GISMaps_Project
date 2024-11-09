/* 
 * Copyright 2024 University of Toronto
 *
 * Permission is hereby granted, to use this software and associated 
 * documentation files (the "Software") in course work at the University 
 * of Toronto, or for personal use. Other uses are prohibited, in 
 * particular the distribution of the Software either publicly or to third 
 * parties.
 *
 * The above copyright notice and this permission notice shall be included in 
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#include "global.h"
#include <iostream>
#include <cmath>
#include "m1.h"
#include "StreetsDatabaseAPI.h"
#include "OSMDatabaseAPI.h"
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <cctype>
#include "ezgl/application.hpp"
#include "ezgl/graphics.hpp"


// loadMap will be called with the name of the file that stores the "layer-2"
// map data accessed through StreetsDatabaseAPI: the street and intersection 
// data that is higher-level than the raw OSM data). 
// This file name will always end in ".streets.bin" and you 
// can call loadStreetsDatabaseBIN with this filename to initialize the
// layer 2 (StreetsDatabase) API.
// If you need data from the lower level, layer 1, API that provides raw OSM
// data (nodes, ways, etc.) you will also need to initialize the layer 1 
// OSMDatabaseAPI by calling loadOSMDatabaseBIN. That function needs the 
// name of the ".osm.bin" file that matches your map -- just change 
// ".streets" to ".osm" in the map_streets_database_filename to get the proper
// name.
double cosineLaw(LatLon A, LatLon B, LatLon C);

std::string getOSMWayValue(OSMID osm_id, std::string key);

extern structIntersection_data Intersection_data;

std::vector<structIntersection_data> intersections; 
std::unordered_map<std::string, std::vector<int>> partialStreetNameStreets;
std::unordered_map<std::string, StreetIdx> streetNamesAndIDs;
std::vector<double> streetLengths; //global variable for every street and its corresponding length
std::vector<std::vector<StreetSegmentIdx>> intersectionStreetSegments; //global variable for finding street segments of an intersection
std::vector<LatLon> intersectionLatLons; //global variable for every intersection's LatLon
std::unordered_map<OSMID,const OSMWay*> OSMWays; //global variable for every OSMWay* attached to its OSMID
std::unordered_map<OSMID,const OSMNode*> OSMNodes;//global variable for every OSMNode* attached to its OSMID
std::unordered_map<StreetIdx,std::vector<IntersectionIdx>> streetIntersections;//global variable of intersections for each street 
std::vector<double> streetSegmentTravelTime; //global variable to store the travel time of every street segment
std::unordered_map<OSMID, std::vector<std::pair<std::string, std::string>>> OSMNodeTag; //global variable
std::unordered_map<OSMID, std::vector<std::pair<std::string, std::string>>> OSMWayTag; 
std::vector<structPOI_data> POI_data; //holds data for all POI

//std::vector <Node> nodes; // store info about nodes; intersection_ID; vector of some data type Node, 1 for each intersection in graph(from src to dest)


//different categories for street segments 
std::vector<StreetIdx> primarySegments;
std::vector<StreetIdx> motorwaySegments;
std::vector<StreetIdx> linkSegments;
std::vector<StreetIdx> secondarySegments;
std::vector<StreetIdx> otherSegments; 
 
//different categories for features 
std::vector<structFeature_data> closedPolygon;
std::vector<structFeature_data> openPolygon; 
std::vector<structFeature_data> singlePoint;

std::vector<structStreetSegment_data> StreetSegment_data;//holds data for each street segment

double min_lon, max_lon, min_lat, max_lat;
double lat_avrg;
float highest_speedLim;

ezgl::renderer *g;
ezgl::surface *restaurantIcon;
ezgl::surface *educationIcon;
ezgl::surface *destinationIcon;
ezgl::surface *dotIcon; 

    //max coordinates for intersections

bool loadMap(std::string map_streets_database_filename) {
  std::cout << "loadMap: " << map_streets_database_filename << std::endl;
  bool load_successful = false;
  if(map_streets_database_filename.length()<13){ //filename too short
    return false;
  }
  if(map_streets_database_filename.substr(map_streets_database_filename.length()-12,map_streets_database_filename.length())!= ".streets.bin"){ //last 12 characters are not .streets.bin
    return false;
  }
    load_successful = loadStreetsDatabaseBIN(map_streets_database_filename); //load streetsdatabase
    if(load_successful){
      std::string OSMFilename = map_streets_database_filename.substr(0, map_streets_database_filename.length()-12); //take out .streets.bin 
      OSMFilename += ".osm.bin"; //replace with .osm.bin
      load_successful = loadOSMDatabaseBIN(OSMFilename); //load osmdatabase
    }else{
      return false;
    }
    //resizing and reserving memory for global variables
    intersectionStreetSegments.resize(getNumIntersections()); 
    intersectionLatLons.resize(getNumIntersections());
    streetSegmentTravelTime.reserve(getNumStreetSegments());
    OSMWays.reserve(getNumberOfWays());
    OSMNodes.reserve(getNumberOfNodes()); 
    streetLengths.resize(getNumStreets());
    intersections.resize(getNumIntersections());
    POI_data.resize(getNumPointsOfInterest());

    //nodes.resize(getNumIntersections()); //resizing and reserving memory for nodes

    //features.resize(getNumFeatures());
    singlePoint.resize(getNumFeatures()); 
    closedPolygon.resize(getNumFeatures()); 
    openPolygon.resize(getNumFeatures());  

    StreetSegment_data.resize(getNumStreetSegments());

    max_lon = getIntersectionPosition(0).longitude();
    min_lon = max_lon;
    max_lat = getIntersectionPosition(0).latitude();
    min_lat = max_lat;

//load icons
    restaurantIcon = g->load_png("libstreetmap/resources/restaurant.png");
    educationIcon = g->load_png("libstreetmap/resources/medicine.png");
    destinationIcon = g->load_png("libstreetmap/resources/location.png");
    dotIcon = g->load_png("libstreetmap/resources/dot.png");

     for(int intersection = 0;intersection < getNumIntersections(); ++intersection){ //cycle all intersections
      intersectionLatLons[intersection] = getIntersectionPosition(intersection);
      intersections[intersection].name = getIntersectionName(intersection);
      //intersections[intersection].xy_loc.x = x_from_lon(getIntersectionPosition(intersection).longitude()); //save its x,y location to the intersections struct
      //intersections[intersection].xy_loc.y = y_from_lat(getIntersectionPosition(intersection).latitude());

      //nodes[intersection].position = getIntersectionPosition(intersection);

        max_lat = std::max(max_lat, getIntersectionPosition(intersection).latitude()); //get the maximum and minimum lat and lon of every intersection
        min_lat = std::min(min_lat, getIntersectionPosition(intersection).latitude());
        max_lon = std::max(max_lon, getIntersectionPosition(intersection).longitude());
        min_lon = std::min(min_lon, getIntersectionPosition(intersection).longitude());
        lat_avrg += kDegreeToRadian*getIntersectionPosition(intersection).latitude();
      
      for(int i = 0; i < getNumIntersectionStreetSegment(intersection);++i){ //cycle all street segments of intersection
        int ss_id = getIntersectionStreetSegment(i, intersection);
        intersectionStreetSegments[intersection].push_back(ss_id); //add street segment to its intersection
        StreetSegmentInfo ss_info = getStreetSegmentInfo(ss_id);
        streetIntersections[ss_info.streetID].push_back(intersection); //add intersection to its street
      }
    }
    lat_avrg /= getNumIntersections();
    for(int intersection = 0;intersection < getNumIntersections(); ++intersection){ //cycle all intersections
      intersections[intersection].xy_loc.x = x_from_lon(getIntersectionPosition(intersection).longitude()); //save its x,y location to the intersections struct
      intersections[intersection].xy_loc.y = y_from_lat(getIntersectionPosition(intersection).latitude());
    }

    for(int way = 0; way < getNumberOfWays(); ++way){ //initialize OSMWays
      const OSMWay* OSMWay = getWayByIndex(way); //get every way's OSMWay*
     OSMID wayID = OSMWay->id(); //find every OSMWay*'s OSMID
     OSMWays[wayID] = OSMWay; //attach OSMWay* to OSMID
     std::vector<std::pair<std::string, std::string >> OSMTag;
      int tags = getTagCount(OSMWay); //num of tags for node 
      OSMTag.resize(tags);

      for (int tag = 0; tag < tags; tag++){
       OSMTag[tag] = getTagPair(OSMWay, tag); //get tagpair for every tag
      }

      OSMWayTag[wayID] = OSMTag; //attach tag pair to the OSMID
    }

    for(int node = 0; node < getNumberOfNodes(); ++node){
      const OSMNode* OSMNode = getNodeByIndex(node); //get every node's OSMNode*
     OSMID nodeID = OSMNode->id(); //find every OSMNode*'s OSMID
     OSMNodes[nodeID] = OSMNode; //attach OSMNode* to OSMID
     std::vector<std::pair<std::string, std::string >> OSMTag;
      int tags = getTagCount(OSMNode); //num of tags for node 
      OSMTag.resize(tags);

      for (int tag = 0; tag < tags; tag++){
       OSMTag[tag] = getTagPair(OSMNode, tag); //get tagpair for every tag
      }

      OSMNodeTag[nodeID] = OSMTag; //attach tag pair to the OSMID
    }

    for(int streetSegment = 0; streetSegment < getNumStreetSegments();streetSegment++){ //cycle all street segments
      StreetSegmentInfo streetSegInfo = getStreetSegmentInfo(streetSegment);
      double travelTime = findStreetSegmentLength(streetSegment)/streetSegInfo.speedLimit; //find traveltime
      streetSegmentTravelTime.push_back(travelTime); //add traveltime to vector
      streetLengths[streetSegInfo.streetID]+= findStreetSegmentLength(streetSegment);
    }

    for(StreetIdx street = 0; street < getNumStreets(); ++street){
      std::string name = getStreetName(street);
      streetNamesAndIDs[name] = street;
      name.erase(std::remove_if(name.begin(),name.end(), ::isspace),name.end()); //erase spaces
      std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {return std::tolower(c);}); //change to lowercase
      int length = getStreetName(street).size();
      for(int i = 1; i < length+1; ++i){
      std::string temporary = name.substr(0, i);
      partialStreetNameStreets[temporary].push_back(street);
      }
    }
    
    highest_speedLim = 0;
    
    for(StreetIdx street = 0; street < getNumStreetSegments(); ++street){ //loop through all street segments
      StreetSegmentInfo ss_info = getStreetSegmentInfo(street); // store all information of each SS into a struct

      OSMID osm_id = ss_info.wayOSMID; 
      LatLon from = getIntersectionPosition(ss_info.from);
      LatLon to = getIntersectionPosition(ss_info.to);
      bool oneWay = ss_info.oneWay;
      int numCurvePoints = ss_info.numCurvePoints; 
      float speedLim = ss_info.speedLimit;
      int StreetID = ss_info.streetID; // index of street this belongs to
      double ss_length = findStreetSegmentLength(street);

      if(speedLim > highest_speedLim){
        highest_speedLim = speedLim; 
      }

      StreetSegment_data[street].ID = StreetID; // stores corressponding street ID for each SS
      StreetSegment_data[street].name = getStreetName(StreetID); // get the name from each ID
      StreetSegment_data[street].oneWay = oneWay;
      StreetSegment_data[street].length = ss_length; 
      StreetSegment_data[street].speedLimit = speedLim; 

      StreetSegment_data[street].numCurvePoints = numCurvePoints; 
   
         if(numCurvePoints == 0){ //segment goes straight from one intersection to next
        StreetSegment_data[street].xy_points.resize(2);
        StreetSegment_data[street].xy_points[0].x = x_from_lon(from.longitude()); //save xy coordinate of start intersection
        StreetSegment_data[street].xy_points[0].y = y_from_lat(from.latitude());
        StreetSegment_data[street].xy_points[1].x = x_from_lon(to.longitude()); //save xy coordinate of end intersection
        StreetSegment_data[street].xy_points[1].y = y_from_lat(to.latitude());
      } else{ //segment has curve points 
        StreetSegment_data[street].xy_points.resize(numCurvePoints + 2); //+2 is including the start point and end point 
        StreetSegment_data[street].xy_points[0].x = x_from_lon(from.longitude()); //save xy coordinate of start intersection
        StreetSegment_data[street].xy_points[0].y = y_from_lat(from.latitude());
        for(int i = 1; i < numCurvePoints + 1; i++){
          StreetSegment_data[street].xy_points[i].x = x_from_lon(getStreetSegmentCurvePoint(i - 1, street).longitude()); //save xy coordinate of curve points
          StreetSegment_data[street].xy_points[i].y = y_from_lat(getStreetSegmentCurvePoint(i - 1, street).latitude());
        }
        StreetSegment_data[street].xy_points[numCurvePoints + 1].x = x_from_lon(to.longitude()); //save xy coordinate of end intersection
        StreetSegment_data[street].xy_points[numCurvePoints + 1].y = y_from_lat(to.latitude());
      }

      if(getOSMWayValue(osm_id, "highway") == "primary"){ //find what type of road it is from OSM
        StreetSegment_data[street].type = "primary"; //categorize different types with different colours and widths
        StreetSegment_data[street].width = 7;
        StreetSegment_data[street].colour = ezgl::color(254, 128, 129); //red
        StreetSegment_data[street].colour_darkMode = ezgl::color(135, 89, 89); //dark red

      } else if(getOSMWayValue(osm_id, "highway") == "motorway"){
        StreetSegment_data[street].type = "motorway";
        StreetSegment_data[street].width = 7;
        StreetSegment_data[street].colour = ezgl::color(238, 185, 76); //orange
        StreetSegment_data[street].colour_darkMode = ezgl::color(166, 136, 104); //dark orange

      } else if(getOSMWayValue(osm_id, "highway") == "motorway_link"){
        StreetSegment_data[street].type = "motorway_link"; 
        StreetSegment_data[street].width = 5.5;
        StreetSegment_data[street].colour = ezgl::color(238, 185, 76); //orange 
        StreetSegment_data[street].colour_darkMode = ezgl::color(166, 136, 104); //dark orange

      } else if(getOSMWayValue(osm_id, "highway") == "secondary"){
        StreetSegment_data[street].type = "secondary";
        StreetSegment_data[street].width = 5;
        StreetSegment_data[street].colour = ezgl::color(242, 238, 233); //very light grey
        StreetSegment_data[street].colour_darkMode = ezgl::color(0x8C, 0x8C, 0x8C); //dark grey

      } else{
        StreetSegment_data[street].type = "other";
        StreetSegment_data[street].width = 5;
        StreetSegment_data[street].colour = ezgl::color(242, 238, 233); //very light grey
        StreetSegment_data[street].colour_darkMode = ezgl::color(0x8C, 0x8C, 0x8C); //dark grey
        
      } 
    } 

    for(POIIdx poi = 0; poi < getNumPointsOfInterest(); poi++){//get all data for POI
      POI_data[poi].xy_points.x = x_from_lon(getPOIPosition(poi).longitude()); //find xy coordinate of POI
      POI_data[poi].xy_points.y = y_from_lat(getPOIPosition(poi).latitude());
      POI_data[poi].type = getPOIType(poi);
      POI_data[poi].name = getPOIName(poi); 
    }
/*
    for (int i = 0; i < getNumIntersections(); i++) { // store all info of intersections into Nodes struct
      nodes[i].position = getIntersectionPosition(i);
    }*/

     for(FeatureIdx feature = 0; feature < getNumFeatures(); feature++){

        int numPoints = getNumFeaturePoints(feature);
        FeatureType featureType = getFeatureType(feature); 
        std::string featureName = getFeatureName(feature);
        ezgl::color featureColour, featureColour_darkMode;
        double x_total = 0, y_total = 0; //use to calculate the center of a feature 
        double x_coord, y_coord;

        switch (featureType){ //categorize the features into dif colours based on their type
          case 1: //parks
            featureColour = {186, 224, 199}; //light green
            featureColour_darkMode = {24, 87, 76}; //pale green
            break; 
          case 2: //beachs
            featureColour = { 254, 255, 225}; //light yellow
            featureColour_darkMode = {145, 134, 105}; //darker beige
            break; 
          case 3: //lakes
            featureColour = { 166, 199, 232}; //light blue
            featureColour_darkMode = {66, 79, 135}; //dark blue
            break;  
          case 4:  //rivers
            featureColour = { 166, 199, 232}; //light blue, same colour as lake
            featureColour_darkMode = {66, 79, 135}; //dark blue
            break; 
          case 5: //islands
            featureColour = { 0xBF, 0xBF, 0xBF}; //grey, same as backdrop
            featureColour_darkMode = {25, 35, 7}; //dark purple

            break; 
          case 6: //buildings
            featureColour = { 0x8C, 0x8C, 0x8C}; //light grey
            featureColour_darkMode = {64, 64, 94}; //dark purple
            break; 
          case 7: //greenspace
            featureColour = { 159, 213, 176}; //dark green
            featureColour_darkMode = { 13, 71, 62};//darker green
            break; 
          case 8: //golfcourse
            featureColour = { 179, 220, 191}; 
            featureColour_darkMode = {24, 87, 76}; //pale green
            break;
          case 9: //streams
            featureColour = { 166, 199, 232}; //light blue, same colour as lake
            featureColour_darkMode = {66, 79, 135}; //dark blue
            break; 
          case 10: //glacier
            featureColour = {87, 118, 145}; //blue
            featureColour_darkMode = {86, 104, 117}; //darker blue
            break;
          default: //unknown
            featureColour = ezgl::GREY_75;
            break;
        }  
         /*if(numPoints == 1){ //feature is a single point
          xy_coord.x = x_from_lon(getFeaturePoint(0, feature).longitude());
          xy_coord.y = y_from_lat(getFeaturePoint(0, feature).latitude()); 
          singlePoint[feature].xy_loc = xy_coord; 
          singlePoint[feature].colour = featureColour; 
          singlePoint[feature].name = featureName; 
        } */
        if (numPoints > 1) { //check if feature is a closed or open polygon
          if(getFeaturePoint(0, feature) == getFeaturePoint(numPoints - 1, feature)){ //if first coord = last coord, closed polygon
            closedPolygon[feature].name = featureName; 
            closedPolygon[feature].numFeaturePoints = numPoints;
            closedPolygon[feature].colour = featureColour;
            closedPolygon[feature].colour_darkMode = featureColour_darkMode;  
            closedPolygon[feature].xy_points.resize(numPoints);
            closedPolygon[feature].area = findFeatureArea(feature);
            for(int i = 0; i < numPoints; i++){ //add coords for all feature points to vec
              x_coord = x_from_lon(getFeaturePoint(i, feature).longitude());
              y_coord = y_from_lat(getFeaturePoint(i, feature).latitude());
              closedPolygon[feature].xy_points[i].x = x_coord; //get x coordinate 
              closedPolygon[feature].xy_points[i].y = y_coord; //get y coordinate

              //find min and max values to calculate xy coordinates for center of feature 
              x_total = x_coord + x_total; 
              y_total = y_coord + y_total;
            }

            //calculate the values for the center of the features
            closedPolygon[feature].x_center = x_total/numPoints;
            closedPolygon[feature].y_center = y_total/numPoints;

          } else{ //open polygon/line feature
            openPolygon[feature].name = featureName; 
            openPolygon[feature].numFeaturePoints = numPoints;
            openPolygon[feature].colour = featureColour;
            openPolygon[feature].xy_points.resize(numPoints);
            openPolygon[feature].colour_darkMode = featureColour_darkMode; 
            for(int i = 0; i < numPoints; i++){ //add coords for all feature points to vec
              x_coord = x_from_lon(getFeaturePoint(i, feature).longitude());
              y_coord = y_from_lat(getFeaturePoint(i, feature).latitude());
              openPolygon[feature].xy_points[i].x = x_coord; //get x coordinate
              openPolygon[feature].xy_points[i].y = y_coord; //get y coordinate
              
              //find min and max values to calculate xy coordinates for center of feature 
              x_total = x_coord + x_total; 
              y_total = y_coord + y_total;
            }

            //calculate the values for the center of the features
            openPolygon[feature].x_center = x_total/numPoints;
            openPolygon[feature].y_center = y_total/numPoints; 
          }
        } 
      }   
    return load_successful;
    //return true;
}                               

void closeMap() {
  closeStreetDatabase();
  closeOSMDatabase();
  streetLengths.clear();
  intersections.clear();
  partialStreetNameStreets.clear();
  streetNamesAndIDs.clear();
  intersectionStreetSegments.clear();
  intersectionLatLons.clear();
  OSMWays.clear();
  OSMNodes.clear();
  streetIntersections.clear();
  streetSegmentTravelTime.clear();
  OSMNodeTag.clear();
  OSMWayTag.clear();
  POI_data.clear();
  primarySegments.clear();
  motorwaySegments.clear();
  linkSegments.clear();
  secondarySegments.clear();
  otherSegments.clear();
  closedPolygon.clear();
  openPolygon.clear();
  singlePoint.clear();
  StreetSegment_data.clear();

  //nodes.clear();

  min_lon = 0;
  max_lon = 0;
  min_lat = 0;
  max_lat = 0; 
  lat_avrg = 0;
  g->free_surface(restaurantIcon);
  g->free_surface(educationIcon);
  g->free_surface(destinationIcon);
  g->free_surface(dotIcon);
      
      // return;
      
}

float x_from_lon(float lon){
return kEarthRadiusInMeters*lon*kDegreeToRadian*cos(lat_avrg);
}

float y_from_lat(float lat){
  return kEarthRadiusInMeters*lat*kDegreeToRadian;
}

float lon_from_x(float x) {
  return (x/(kEarthRadiusInMeters * kDegreeToRadian * cos(lat_avrg)));
}

float lat_from_y(float y) {
  return (y/(kEarthRadiusInMeters * kDegreeToRadian));
}

double findDistanceBetweenTwoPoints(LatLon point_1, LatLon point_2) { //sissi 
  double distance, lat, lon, lat2, lon2, lat_avg; 
  double x, y, x2, y2; 

  //for coordinates (x, y): x = R*lon*cos(lat_avg), Y = R*lat
  //lat_avg = average latitude of area being mapped, R = radius of the Earth, lat lon in radians 

  //convert lat lon to radians 
  lat = point_1.latitude()*kDegreeToRadian;
  lon = point_1.longitude()*kDegreeToRadian;
  lat2 = point_2.latitude()*kDegreeToRadian;
  lon2 = point_2.longitude()*kDegreeToRadian; 

  //find lat_avg
  lat_avg = (lat + lat2)/2.0;

  //x = R*lon*cos(lat_avg)
  x = kEarthRadiusInMeters*lon*cos(lat_avg);
  x2 = kEarthRadiusInMeters*lon2*cos(lat_avg);

  //y = R*lat 
  y = kEarthRadiusInMeters*lat;
  y2 = kEarthRadiusInMeters*lat2;

  distance = sqrt(pow(y2-y, 2) + pow(x2-x, 2));
  return distance; 
}

// Returns the length of the given street segment in meters.
// Speed Requirement --> moderate
double findStreetSegmentLength(StreetSegmentIdx street_segment_id) { 

  double distance = 0.0;
  LatLon l_id1, l_id2;

  std::vector <LatLon> val; // vector type LatLon (takes in LatLon types), named val; to store all points on a SS (used if there is a curve on street)

  StreetSegmentInfo ss_id = getStreetSegmentInfo(street_segment_id); //access struct info for SS of given ID

  l_id1 = getIntersectionPosition(ss_id.from); //convert intersection A to LatLon
  l_id2 = getIntersectionPosition(ss_id.to); //convert intersection B to LatLon
 
  if (ss_id.numCurvePoints == 0) { // no curve points; only a straight line 
    return findDistanceBetweenTwoPoints(l_id1, l_id2);
  }
  else {
    val.push_back(l_id1); // make the first element in vector = Intersection A
    for (int i = 0; i < ss_id.numCurvePoints; ++i) { // loop through curve points
      LatLon point = getStreetSegmentCurvePoint(i, street_segment_id); // convert the curve points index into LatLon      
      val.push_back(point); // put all LatLon values in vector
    }
    val.push_back(l_id2); // when done with curve, push in Intersection B
    for (int i = 0; i < val.size()-1; ++i) { 
      distance = distance + findDistanceBetweenTwoPoints(val[i], val[i+1]);  // to calculate distance b/w 2 adjacent points
    } // return when you have (size - 1) to avoid reaching end of vector error
    return distance;
  }
}

double findStreetSegmentTravelTime(StreetSegmentIdx street_segment_id){
 return streetSegmentTravelTime[street_segment_id];
}

double findAngleBetweenStreetSegments(StreetSegmentIdx src_street_segment_id, StreetSegmentIdx dst_street_segment_id){
//retrieve street information
  StreetSegmentInfo ss_info_src = getStreetSegmentInfo (src_street_segment_id);
  StreetSegmentInfo ss_info_dst = getStreetSegmentInfo (dst_street_segment_id);
  int numCurvePoints = ss_info_src.numCurvePoints; 
  int numCurvePoints2 = ss_info_dst.numCurvePoints;
 
  LatLon a_to, a_from, b_to, b_from, c;

  //find points a,b needed to apply cosine law
  if(numCurvePoints == 0){ //no curve points in ss1
    a_to = getIntersectionPosition(ss_info_src.from);
    a_from = getIntersectionPosition(ss_info_src.to);
  } else{ //there are curve points in ss1 
    a_from = getStreetSegmentCurvePoint(0, src_street_segment_id);
    a_to = getStreetSegmentCurvePoint(numCurvePoints-1 , src_street_segment_id);
  }

  if(numCurvePoints2 == 0){ //no curve points in ss2
    b_to = getIntersectionPosition(ss_info_dst.from);
    b_from = getIntersectionPosition(ss_info_dst.to);
  } else{ //there are curve points in ss2
    b_from = getStreetSegmentCurvePoint(0, dst_street_segment_id);
    b_to = getStreetSegmentCurvePoint(numCurvePoints2-1 , dst_street_segment_id);
  }
 
  //get cosARG
  if(ss_info_src.from == ss_info_dst.to){
    c = getIntersectionPosition(ss_info_src.from);
    return cosineLaw(a_from, b_to, c);

  } else if(ss_info_src.from == ss_info_dst.from){
    c = getIntersectionPosition(ss_info_src.from);
    return cosineLaw(a_from, b_from, c);

  } else if(ss_info_src.to == ss_info_dst.to){
    c = getIntersectionPosition(ss_info_src.to);
    return cosineLaw(a_to, b_to, c);

  } else if(ss_info_src.to == ss_info_dst.from){
    c = getIntersectionPosition(ss_info_src.to);
    return cosineLaw(a_to, b_from, c);

  } else{ //no shared intersections
    return NO_ANGLE; 
  }
}

//return cosArg
double cosineLaw(LatLon A, LatLon B, LatLon C){
  //find lengths to each side of "triangle"
  double a = findDistanceBetweenTwoPoints(C, B);
  double b = findDistanceBetweenTwoPoints(A, C); 
  double c = findDistanceBetweenTwoPoints(A, B);

  //use cos law 
  double argCos = (pow(a,2)+pow(b,2)-pow(c,2))/(2.0*a*b); 

  //edge cases, check if out of range [-1, 1]
  if (argCos >= 1.0 ||argCos <= -1.0){
    return 0; 
  } else { 
    return (M_PI - acos(argCos)); 
  }
}



// Returns true if the two intersections are directly connected, meaning you can
// legally drive from the first intersection to the second using only one
// streetSegment.
// Speed Requirement --> moderate
bool intersectionsAreDirectlyConnected(std::pair<IntersectionIdx, IntersectionIdx> intersection_ids){
  StreetSegmentIdx found;
  if (intersection_ids.first == intersection_ids.second) { // if the 2 IDs are the same
    return true;
  }
  else {
  // first access all SS of 1st intersection (will receive a vector from calling function):
    std::vector <StreetSegmentIdx> a = findStreetSegmentsOfIntersection(intersection_ids.first);
  // then access all SS of 2nd intersection:
    std::vector <StreetSegmentIdx> b = findStreetSegmentsOfIntersection(intersection_ids.second);
  // then check if they share a SS in common; if not return false:
    std::vector <StreetSegmentIdx> c(a.size() + b.size()); //set the size of c (total) to be the size of the 2 combined
    std::vector <StreetSegmentIdx> :: iterator it;
    it = std::set_intersection (a.begin(), a.end(), b.begin(), b.end(), c.begin()); // iterates and compares vectors a and b, and places same SS in c
    c.resize(it - c.begin()); //cuts down extra space
    if (c.size() == 0) { // no common elements because cut down everything, so empty vector c
      return false;
    }
    else { // found a common SS
      found = c[0]; // set found to be the SS number 
       StreetSegmentInfo ss_id = getStreetSegmentInfo(found); //access struct info for SS of found ID
      if (!ss_id.oneWay) { 
        return true; // return true if 2 way
      }
      else {
        if((intersection_ids.first == ss_id.from) && (intersection_ids.second == ss_id.to)){ // return true if the inputs match from and to 
          return true;
        }
        else {
          return false;
        }
      }
    }
  }
} 

IntersectionIdx findClosestIntersection(LatLon my_position){
  double closestDistance;
  IntersectionIdx closestIntersection = 0;
  LatLon intersectionLatLon;
  intersectionLatLon = intersectionLatLons[0]; //store first intersection LatLon 
  closestDistance = findDistanceBetweenTwoPoints(my_position, intersectionLatLon); //store the distance from first intersection LatLon as closest so far
for (int intersection = 1; intersection < intersectionLatLons.size()-1; ++intersection){ //loop through all intersections
intersectionLatLon= intersectionLatLons[intersection]; 
if(findDistanceBetweenTwoPoints(my_position,intersectionLatLon)<closestDistance){ //if intersection is found to have shorter disance, make that the closest intersection
  closestIntersection = intersection;
  closestDistance = findDistanceBetweenTwoPoints(my_position,intersectionLatLon);
}
}
return closestIntersection;
}

// Returns the street segments that connect to the given intersection.
// Speed Requirement --> high
std::vector<StreetSegmentIdx> findStreetSegmentsOfIntersection(IntersectionIdx intersection_id){
  return intersectionStreetSegments[intersection_id];
}

std::vector<IntersectionIdx> findIntersectionsOfStreet(StreetIdx street_id) { //given in lecture 
 std::sort(streetIntersections[street_id].begin(), streetIntersections[street_id].end());
  auto duplicates = unique(streetIntersections[street_id].begin(),streetIntersections[street_id].end()); //take all duplicates of intersections vector
  streetIntersections[street_id].erase(duplicates, streetIntersections[street_id].end()); //erase from vector
  return streetIntersections[street_id]; //return vector of unique intersections
}

std::vector<IntersectionIdx> findIntersectionsOfTwoStreets(std::pair<StreetIdx, StreetIdx> street_ids){
  StreetIdx street = street_ids.first; 
  StreetIdx street2 = street_ids.second; 
  std::vector<IntersectionIdx> inter_vec; 
  std::vector<IntersectionIdx> street_inter = findIntersectionsOfStreet(street); 
  std::vector<IntersectionIdx> street2_inter = findIntersectionsOfStreet(street2);

  for(int i = 0; i < street_inter.size(); ++i){
    for(int a = 0; a < street2_inter.size(); ++a){
      if(street_inter[i] == street2_inter[a]){
        inter_vec.push_back(street2_inter[a]); 
        break; 
      }
    }
  }
  return inter_vec; 
}

// Returns all street ids corresponding to street names that start with the
// given prefix.
// The function should be case-insensitive to the street prefix.
// The function should ignore spaces.
//  For example, both "bloor " and "BloOrst" are prefixes to 
//  "Bloor Street East".

// If no street names match the given prefix, this routine returns an empty
// (length 0) vector.
// You can choose what to return if the street prefix passed in is an empty
// (length 0) string, but your program must not crash if street_prefix is a
// length 0 string.

// Speed Requirement --> high
std::vector<StreetIdx> findStreetIdsFromPartialStreetName(std::string street_prefix){ 
 street_prefix.erase(std::remove_if(street_prefix.begin(),street_prefix.end(), ::isspace),
 street_prefix.end()); //remove spaces
 std::transform(street_prefix.begin(), street_prefix.end(), street_prefix.begin(), [](unsigned char c) {return std::tolower(c);}); //change to lowercase
 std::sort(partialStreetNameStreets[street_prefix].begin(), partialStreetNameStreets[street_prefix].end());
  auto duplicates = unique(partialStreetNameStreets[street_prefix].begin(),partialStreetNameStreets[street_prefix].end()); //take all duplicates of intersections vector
  partialStreetNameStreets[street_prefix].erase(duplicates, partialStreetNameStreets[street_prefix].end()); //erase from vector
 return partialStreetNameStreets[street_prefix];
}

double findStreetLength(StreetIdx street_id){
  return streetLengths[street_id];
}

POIIdx findClosestPOI(LatLon my_position, std::string poi_name) {
  std::string compare_name;
  double shortest_distance = 100000000;
  double distance;
  LatLon poi_position; 
  POIIdx closest_POI = -1;

   for(POIIdx i = 0; i < getNumPointsOfInterest(); i++){
    compare_name = getPOIName(i);
    if(compare_name == poi_name){
      poi_position = getPOIPosition(i); 
      distance = findDistanceBetweenTwoPoints(my_position, poi_position);
      if(distance < shortest_distance){ 
        closest_POI = i;
        shortest_distance = distance; 
      }
    }
  }
  return closest_POI; 
}

// Returns the area of the given closed feature in square meters.
// Assume a non self-intersecting polygon (i.e. no holes).
// Return 0 if this feature is not a closed polygon.
// Speed Requirement --> moderate
double findFeatureArea(FeatureIdx feature_id) {

  std:: vector <double> x_p;
  std:: vector <double> y_p;
  std:: vector <LatLon> p;

  double lat, lon, lat_avg, x, y, area;
  double total_lat = 0.0; 
  double total_area = 0.0;

  int points = getNumFeaturePoints(feature_id);

  for(int i = 0; i < points; ++i) {
    LatLon convert = getFeaturePoint(i, feature_id); // get all points of feature into LatLon value
    p.push_back(convert); // store all LatLon points in the vector
  }
  if (p[0] == (p[points-1])) {  // if first location of point in order, == last location of point ; closedPolygon
    for (int i = 0; i < points; ++ i) {
      lat = p[i].latitude()*kDegreeToRadian; // i.e. take point 0 in vector p, take its latitude and convert to rads,and place that value into lat. //convert lat to radians
      total_lat = total_lat + lat;      //  find the LatAvg of all feature points
      x_p.push_back(lat); // store lats into one vector in order
      
      lon = p[i].longitude()*kDegreeToRadian;  //convert lon to radians
      y_p.push_back(lon); // store lons into another vector in order
    }

    lat_avg = total_lat / points;   // find lat_avg of all feature points:

    for (int i = 0; i < points; ++ i) {
      x = kEarthRadiusInMeters*(y_p[i])*cos(lat_avg); //x = R*lon*cos(lat_avg) ; lon = y_p[i]
      y = kEarthRadiusInMeters*(x_p[i]); //y = R*lat  ; lat = x_p[i]

      x_p[i] = x; // replace lat value in x_p[i] with new calculated x-cor val
      y_p[i] = y; //replace lon value in y_p[i] with new calculated y-cor val
    }
    // now, x_p and y_p has (x,y) coord of each point accordingly; index(point) = x_p(point) , y_p(point)
    // in a loop, calculate area between consecutive points:
    for(int i = 0; i < (points - 1); ++i) {
      area = ((y_p[i+1] - y_p[i]) * ((x_p[i+1]+ x_p[i])/2));
      total_area = total_area + area;
    }
    if (total_area < 0.0) {
      return fabs(total_area);
    }
    return total_area;
  }
  else {   // if first location of point in order, != last location of point ; polyline, return 0
    return 0.0;
  }
}

double findWayLength(OSMID way_id) { //uses OSMDatabase.h
    double length = 0;
    const std::vector<OSMID>& Nodes = getWayMembers(OSMWays[way_id]); //all the nodes of the way's OSMIDs in a vector
    for(int i = 0; i< Nodes.size()-1; ++i){ //for every node
      LatLon node1LatLon = getNodeCoords(OSMNodes[Nodes[i]]); //store latlon of current node
      LatLon node2LatLon = getNodeCoords(OSMNodes[Nodes[i+1]]); //store latlon of next node
      length += findDistanceBetweenTwoPoints(node1LatLon, node2LatLon); //calculate distance and add to total
    }
    return length;
}

std::string getOSMNodeTagValue(OSMID osm_id, std::string key) {  // uses OSMDatabase.h
  auto node = OSMNodeTag.find(osm_id);
  if(node != OSMNodeTag.end()){
    for(const auto& pair: node->second){
      if(pair.first == key){
        return pair.second; 
      }
    }
  }
  return ""; 
}

std::string getOSMWayValue(OSMID osm_id, std::string key) {  // uses OSMDatabase.h
  auto way = OSMWayTag.find(osm_id);
  if(way != OSMWayTag.end()){
    for(const auto& pair: way->second){
      if(pair.first == key){
        return pair.second; 
      }
    }
  }
  return ""; 
}
