CXXFLAGS = -Wall -Wextra -Wno-long-long -Wno-unused-result -Wno-sign-compare -Wno-delete-non-virtual-dtor -Wno-terminate -std=c++17 -fPIC -O3

INCL = -I. -I.. -I../include -I../libslic3r -I/usr/include/eigen3 -I/usr/local/include
LIBS = -lz3
TBB_LIBS = -ltbb -ltbbmalloc
EXTRA_LIBS = -lgecodedriver -lgecodesearch -lgecodeminimodel -lgecodeint -lgecodekernel -lgecodesupport -lgecodegist


all:		seq_test_sequential seq_test_polygon seq_test_preprocess seq_test_interface seq_test_arrangement exp_gecode_arrange exp_z3_arrange exp_cegar_seq


prusaparts.o:			../include/seq_defs.hpp ../include/seq_version.hpp prusaparts.hpp prusaparts.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oprusaparts.o prusaparts.cpp

seq_test_arrangement.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp seq_test_arrangement.hpp seq_test_arrangement.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_test_arrangement.o seq_test_arrangement.cpp

exp_gecode_arrange.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp exp_gecode_arrange.hpp exp_gecode_arrange.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oexp_gecode_arrange.o exp_gecode_arrange.cpp

seq_test_polygon.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp seq_test_polygon.hpp seq_test_polygon.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_test_polygon.o seq_test_polygon.cpp

seq_test_preprocess.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp seq_test_preprocess.hpp seq_test_preprocess.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_test_preprocess.o seq_test_preprocess.cpp

exp_cegar_seq.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp exp_cegar_seq.hpp exp_cegar_seq.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oexp_cegar_seq.o exp_cegar_seq.cpp

seq_test_interface.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_interface.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp seq_test_interface.hpp seq_test_interface.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_test_interface.o seq_test_interface.cpp

seq_test_sequential.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp seq_test_sequential.hpp seq_test_sequential.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oseq_test_sequential.o seq_test_sequential.cpp

exp_z3_arrange.o:		../include/seq_defs.hpp ../include/seq_version.hpp ../include/seq_utilities.hpp ../include/seq_preprocess.hpp ../include/seq_sequential.hpp ../include/seq_interface.hpp exp_z3_arrange.hpp exp_z3_arrange.cpp
				g++ $(CXXFLAGS) $(INCL) -c -oexp_z3_arrange.o exp_z3_arrange.cpp

seq_test_slicer:		seq_test_slicer.o
				g++ $(CXXFLAGS) $(INCL)

seq_test_arrangement:		seq_test_arrangement.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) -L/usr/local/lib -L/usr/lib seq_test_arrangement.o -oseq_test_arrangement  $(LIBS) $(EXTRA_LIBS)

exp_gecode_arrange:		exp_gecode_arrange.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) -L/usr/local/lib -L/usr/lib exp_gecode_arrange.o -oexp_gecode_arrange  $(LIBS) $(EXTRA_LIBS)

seq_test_polygon:		seq_test_polygon.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) -L/usr/local/lib seq_test_polygon.o ../src/seq_preprocess.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o ../src/seq_sequential.o -oseq_test_polygon -L/usr/lib -L../libslic3r -lSlic3r_tiny $(LIBS) $(TBB_LIBS)

seq_test_sequential:		seq_test_sequential.o ../src/seq_sequential.o ../src/seq_preprocess.o prusaparts.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) seq_test_sequential.o ../src/seq_preprocess.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o -oseq_test_sequential -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

exp_z3_arrange:			exp_z3_arrange.o ../src/seq_sequential.o ../src/seq_preprocess.o prusaparts.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) exp_z3_arrange.o ../src/seq_preprocess.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o -oexp_z3_arrange -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)


seq_test_preprocess:		seq_test_preprocess.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) seq_test_preprocess.o ../src/seq_preprocess.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o -oseq_test_preprocess -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

exp_cegar_seq:			exp_cegar_seq.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) exp_cegar_seq.o ../src/seq_preprocess.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o -oexp_cegar_seq -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

seq_test_interface:		seq_test_interface.o ../src/seq_preprocess.o prusaparts.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o
				g++ $(CXXFLAGS) $(INCL) seq_test_interface.o ../src/seq_preprocess.o ../src/seq_sequential.o ../src/seq_interface.o ../src/seq_utilities.o prusaparts.o -oseq_test_interface -L../libslic3r -lSlic3r_tiny -L/usr/local/lib -L/usr/lib $(LIBS) $(TBB_LIBS)

clean:
				rm -f sequential_prusa sequential_decimator seq_test_arrangement exp_gecode_arrange seq_test_sequential exp_z3_arrange seq_test_polygon seq_test_preprocess exp_cegar_seq seq_test_interface *.o *.o *~ Geometry/*.o libslic3r/*.o *~ Geometry/*.a libslic3r/*.a clipper/*.o
