CXXFLAGS = -Wall -Wextra -pedantic -Wno-long-long -Wno-unused-result -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-terminate -std=c++17 -fPIC -O3

INCL = -I. -I./.. -I/usr/include/eigen3 -I/usr/local/include
LIBS = -lgecodedriver -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodekernel -lgecodesupport -lgecodegist -lz3

HEADERS = Polygon.hpp BoundingBox.hpp Line.hpp Geometry.hpp ExPolygon.hpp MultiPoint.hpp Point.hpp Polyline.hpp SVG.hpp Surface.hpp clipper.hpp ClipperUtils.hpp TriangleMesh.hpp ConvexHull.hpp ./Geometry/MedialAxis.hpp ./Geometry/VoronoiOffset.hpp LocalesUtils.hpp ./clipper/clipper.hpp

debug:		all
release:	all

all:	libSlic3r_tiny.a

Polygon.o:	$(HEADERS) Polygon.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oPolygon.o Polygon.cpp

ExPolygon.o:	$(HEADERS) ExPolygon.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oExPolygon.o ExPolygon.cpp

BoundingBox.o:	$(HEADERS) BoundingBox.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oBoundingBox.o BoundingBox.cpp

Line.o:		$(HEADERS) Line.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oLine.o Line.cpp

Geometry.o:	$(HEADERS) Geometry.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oGeometry.o Geometry.cpp

./Geometry/MedialAxis.o:	$(HEADERS) ./Geometry/MedialAxis.cpp
				g++ $(CXXFLAGS) $(INCL) -c -o./Geometry/MedialAxis.o ./Geometry/MedialAxis.cpp

./Geometry/VoronoiOffset.o:	$(HEADERS) ./Geometry/VoronoiOffset.cpp
				g++ $(CXXFLAGS) $(INCL) -c -o./Geometry/VoronoiOffset.o ./Geometry/VoronoiOffset.cpp

MultiPoint.o:	$(HEADERS) MultiPoint.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oMultiPoint.o MultiPoint.cpp

Point.o:	$(HEADERS) Point.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oPoint.o Point.cpp

Polyline.o:	$(HEADERS) Polyline.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oPolyline.o Polyline.cpp

SVG.o:		$(HEADERS) SVG.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oSVG.o SVG.cpp

Surface.o:	$(HEADERS) Surface.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oSurface.o Surface.cpp

LocalesUtils.o:	$(HEADERS) LocalesUtils.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oLocalesUtils.o LocalesUtils.cpp

ConvexHull.o:	$(HEADERS) ConvexHull.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oConvexHull.o ConvexHull.cpp

clipper.o:	$(HEADERS) clipper.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oclipper.o clipper.cpp

./clipper/clipper.o:	$(HEADERS) ./clipper/clipper.cpp
			g++ $(CXXFLAGS) $(INCL) -c -o./clipper/clipper.o ./clipper/clipper.cpp


ClipperUtils.o:	$(HEADERS) ClipperUtils.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oClipperUtils.o ClipperUtils.cpp

ShortestPath.o:	$(HEADERS) ShortestPath.cpp
		g++ $(CXXFLAGS) $(INCL) -c -oShortestPath.o ShortestPath.cpp

libSlic3r_tiny.a:	Polygon.o BoundingBox.o Line.o Geometry.o ExPolygon.o MultiPoint.o Point.o Polyline.o SVG.o Surface.o clipper.o ClipperUtils.o ShortestPath.o ConvexHull.o ./Geometry/MedialAxis.o ./Geometry/VoronoiOffset.o LocalesUtils.o ./clipper/clipper.o
			ar r libSlic3r_tiny.a Polygon.o BoundingBox.o Line.o Geometry.o ExPolygon.o MultiPoint.o Point.o Polyline.o SVG.o Surface.o clipper.o ClipperUtils.o ShortestPath.o ConvexHull.o ./Geometry/MedialAxis.o ./Geometry/VoronoiOffset.o LocalesUtils.o ./clipper/clipper.o


clean:
	rm -f *.a *.o *~ ./Geometry/*.o ./clipper/*.o
