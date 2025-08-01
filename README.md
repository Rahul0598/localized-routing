# Localized Routing

Interactive web map for routing across a city's walkways and streets, using custom GeoJSON data. This script uses: 

1. A* Routing Algorithm - finds the shortest path between two points on the network of sidewalks and roads. It guarantees the shortest route, as long as the heuristic function does not overestimate the remaining distance.
2. RBush Spatial Index - stores all network nodes for instant nearest-neighbor lookup. Finding the nearest network node to the user's click must be fast, as the network may consist of thousands of points. RBush performs sublinear search time even as data size grows, unlike simple linear scans which become slow for large networks. Enables 'snapping' start/end points to the closest node within an appropriate distance, ensuring the route starts/ends logically on the mapped network.

![Screenshot of Routing Map](./demo.png)

