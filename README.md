# Line-Segment-Intersections
Computing the intersection points and convex hull (using Chan's Algorithm) resulting from those intersection points, from given n - line segments in a 2-D plane in C++

### Input 
Following should be the content of the input file
  - First Line denotes the No. of Line Segments, say T
  - Next T lines denote the end-points, space separated x and y coordinates respectively with the top endpoint           entered first and then the lower one
  - Sample Input:
    3
    0 2 1 0
    0 1 3 0
    2 1 0 -2
  - Here Line Segment #1 is Line Segment joining (0,2) and (1,0) so on and so forth.

###Output
Space separated intersection points and on the next line convex hull set
  - Sample output for above sample input:
    (0.6,0.8) (1.63636,0.454545)  <--- Intersection Points
    (1.63636,0.454545) (0.6,0.8) (1.63636,0.454545)  <--- Convex Hull Points

###Compile Instructions
While compiling and executing, please use redirection operator '<' (for input redirection from file) and '>' (for output redirection to file). For example,
  - $g++ lineSegmentIntersection.cpp
  - $./a.out < input > output



