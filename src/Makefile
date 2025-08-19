CXXFLAGS = -Wall -Wextra -Wno-long-long -Wno-unused-result -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-terminate -std=c++17 -fPIC -O3

INCL = -I. -I.. -I./include -I/usr/include/eigen3 -I/usr/local/include -I../include -I../libslic3r
LIBS = -lz3
TBB_LIBS = -ltbb -ltbbmalloc
EXTRA_LIBS = -lgecodedriver -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodekernel -lgecodesupport -lgecodegist

all:		sequential_decimator sequential_prusa


seq_sequential.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_sequential.hpp ../include/seq_preprocess.hpp ../include/seq_interface.hpp seq_sequential.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_sequential.o seq_sequential.cpp

seq_preprocess.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_sequential.hpp ../include/seq_preprocess.hpp ../include/seq_interface.hpp seq_preprocess.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_preprocess.o seq_preprocess.cpp

seq_utilities.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_sequential.hpp ../include/seq_utilities.hpp ../include/seq_interface.hpp seq_utilities.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_utilities.o seq_utilities.cpp

seq_interface.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_sequential.hpp ../include/seq_preprocess.hpp ../include/seq_interface.hpp seq_interface.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_interface.o seq_interface.cpp

sequential_prusa.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp sequential_prusa.hpp sequential_prusa.cpp
				g++ $(CXXFLAGS) $(INCL) -c -osequential_prusa.o sequential_prusa.cpp

sequential_decimator.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp sequential_decimator.hpp sequential_decimator.cpp
				g++ $(CXXFLAGS) $(INCL) -c -osequential_decimator.o sequential_decimator.cpp

sequential_decimator:		sequential_decimator.o seq_preprocess.o seq_sequential.o seq_interface.o seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) sequential_decimator.o seq_preprocess.o seq_sequential.o seq_interface.o seq_utilities.o -osequential_decimator -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

sequential_prusa:		sequential_prusa.o seq_preprocess.o seq_sequential.o seq_interface.o seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) sequential_prusa.o seq_preprocess.o seq_sequential.o seq_interface.o seq_utilities.o -osequential_prusa -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

clean:
				rm -f sequential_prusa sequential_decimator seq_test_arrangement seq_test_sequential seq_test_polygon seq_test_preprocess seq_test_interface *.o *.o *~ Geometry/*.o libslic3r/*.o *~ Geometry/*.a libslic3r/*.a clipper/*.o
