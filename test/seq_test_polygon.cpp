/*================================================================*/
/*
 * Author:  Pavel Surynek, 2023 - 2024
 * Company: Prusa Research
 *
 * File:    seq_test_polygon.cpp
 *
 * Basic polygon tests.
 */
/*================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/ConvexHull.hpp"
#include "libslic3r/SVG.hpp"

#include <z3++.h>

#include "prusaparts.hpp"

#include "seq_defs.hpp"

#include "seq_sequential.hpp"
#include "seq_preprocess.hpp"

#include "seq_test_polygon.hpp"


/*----------------------------------------------------------------*/

using namespace Slic3r;
using namespace Slic3r::Geometry;

using namespace z3;

using namespace Sequential;

#define SCALE_FACTOR 100000


/*----------------------------------------------------------------*/


void test_polygon_1(void)
{ 
    printf("Testing polygon 1 ...\n");

    Polygon polygon_1 = {{-1000000, -1000000}, {1000000, -1000000}, {1000000, 1000000}, {-1000000, 1000000} };

    for (unsigned int i = 0; i < polygon_1.size(); ++i)
    {
	Point point = polygon_1[i];
	printf("%d,%d\n", point.x(), point.y());
    }
    
    printf("Testing polygon 1 ... finished\n");
}


void test_polygon_2(void)
{ 
    printf("Testing polygon 2 ...\n");

    for (unsigned int k = 0; k < PRUSA_PART_POLYGONS.size(); ++k)
    {
	printf("k = %d\n", k);
	
	const Polygon &polygon_1 = PRUSA_PART_POLYGONS[k];
	Polygon hull_1 = convex_hull(polygon_1);
	
	for (unsigned int i = 0; i < polygon_1.size(); ++i)
	{
	    const Point &point = polygon_1[i];
	    printf("poly %d: %d,%d\n", i, point.x(), point.y());
	}
	printf("\n");
	
	for (unsigned int i = 0; i < hull_1.size(); ++i)
	{
	    const Point &point = hull_1[i];
	    printf("hull %d: %d,%d\n", i, point.x(), point.y());
	}

	if (hull_1.size() >= 2)
	{
	    const Point &point_1 = hull_1[0];
	    const Point &point_2 = hull_1[1];

	    Point v = (point_2 - point_1); //.normalized();
	    printf("v: %d,%d\n", v.x(), v.y());
	    cout << v << endl;
	    
	    Point u = v.normalized();	    
	    printf("u: %d,%d\n", u.x(), u.y());
	    cout << u << endl;

	    printf("Ortho:\n");
	    Point n(v.y(), -v.x());
	    cout << n << endl;

	    coord_t d = n.x() * point_1.x() + n.y() * point_1.y();
	    printf("%d\n", d);
	    cout << d << endl;

	    auto is_inside=[&](const Point &p)
	    {
		coord_t d1 = n.x() * p.x() + n.y() * p.y() - d;
		printf("d1: %d\n", d1);

		if (d1 >= 0)
		{
		    return true;
		}
		else
		{
		    return false;
		}
	    };
	    
	    bool ins1 = is_inside(point_1);
	    printf("%s\n", ins1 ? "yes" : "no");
	    bool ins2 = is_inside(point_2);
	    printf("%s\n", ins2 ? "yes" : "no");	    
	    bool ins3 = is_inside(point_1 + point_2);
	    printf("%s\n", ins3 ? "yes" : "no");	    
	    bool ins4 = is_inside(point_1 - point_2);
	    printf("%s\n", ins4 ? "yes" : "no");	    
	}	
	
	getchar();
    }
    
    printf("Testing polygon 2 ... finished\n");
}


int line_count = 4;
Line lines[] = {{Point(100,100), Point(200,200)}, {Point(200,100), Point(100,200)}, {Point(0,0), Point(100,10)}, {Point(50,0), Point(60,100)} };

void test_polygon_3(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 3 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T_parameters(z_context);    
    
    for (int i = 0; i < line_count; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < line_count; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < line_count; ++i)
    {
	string name = "t_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T_parameters.push_back(expr(z_context.real_const(name.c_str())));
    }    
    
    z3::solver z_solver(z_context);

    introduce_LineNonIntersection_explicit(z_solver,
					   z_context,
					   X_positions[0],
					   Y_positions[0],
					   T_parameters[0],
					   lines[0],
					   X_positions[1],
					   Y_positions[1],
					   T_parameters[1],
					   lines[1]);

    introduce_LineNonIntersection_explicit(z_solver,
					   z_context,
					   X_positions[2],
					   Y_positions[2],
					   T_parameters[2],
					   lines[2],
					   X_positions[3],
					   Y_positions[3],
					   T_parameters[3],
					   lines[3]);    

    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    switch (z_solver.check())
    {
    case z3::sat:
    {
	printf("  SATISFIABLE\n");
	break;
    }
    case z3::unsat:	
    {
	printf("  UNSATISFIABLE\n");
	return;
	break;
    }
    case z3::unknown:
    {
	printf("  UNKNOWN\n");
	break;
    }
    default:
    {
	break;
    }
    }   

    z3::model z_model(z_solver.get_model());
    printf("Printing model:\n");
    cout << z_model << "\n";

    finish = clock();    

    printf("Printing interpretation:\n");    
    for (unsigned int i = 0; i < z_model.size(); ++i)
    {
	printf("Variable:%s  ", z_model[i].name().str().c_str());
	
	cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
	double value = z_model.get_const_interp(z_model[i]).as_double();
	
	printf("value: %.3f\n", value);
	
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 3 ... finished\n");    
}


void test_polygon_4(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 4 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T_parameters(z_context);    
    
    for (int i = 0; i < line_count; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < line_count; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < line_count; ++i)
    {
	string name = "t_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T_parameters.push_back(expr(z_context.real_const(name.c_str())));
    }    
    
    z3::solver z_solver(z_context);

    introduce_LineNonIntersection_implicit(z_solver,
					   z_context,
					   X_positions[0],
					   Y_positions[0],
					   T_parameters[0],
					   lines[0],
					   X_positions[1],
					   Y_positions[1],
					   T_parameters[1],
					   lines[1]);

    introduce_LineNonIntersection_implicit(z_solver,
					   z_context,
					   X_positions[2],
					   Y_positions[2],
					   T_parameters[2],
					   lines[2],
					   X_positions[3],
					   Y_positions[3],
					   T_parameters[3],
					   lines[3]);    

    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    switch (z_solver.check())
    {
    case z3::sat:
    {
	printf("  SATISFIABLE\n");
	break;
    }
    case z3::unsat:	
    {
	printf("  UNSATISFIABLE\n");
	return;
	break;
    }
    case z3::unknown:
    {
	printf("  UNKNOWN\n");
	break;
    }
    default:
    {
	break;
    }
    }   

    z3::model z_model(z_solver.get_model());
    printf("Printing model:\n");
    cout << z_model << "\n";

    finish = clock();    

    printf("Printing interpretation:\n");    
    for (unsigned int i = 0; i < z_model.size(); ++i)
    {
	printf("Variable:%s  ", z_model[i].name().str().c_str());
	
	cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
	double value = z_model.get_const_interp(z_model[i]).as_double();
	
	printf("value: %.3f\n", value);
	
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 4 ... finished\n");    
}


int poly_line_count = 4;
Line poly_lines[] = {{Point(100,100), Point(200,100)}, {Point(200,100), Point(200,200)}, {Point(200,200), Point(100,200)}, {Point(100,200), Point(100,100)} };


void test_polygon_5(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 5 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    
    for (int i = 0; i < poly_line_count; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < poly_line_count; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }
    
    z3::solver z_solver(z_context);

    introduce_PointInsideHalfPlane(z_solver,
				   X_positions[0],
				   Y_positions[0],
				   X_positions[1],
				   Y_positions[1],
				   poly_lines[0]);

    introduce_PointInsideHalfPlane(z_solver,
				   X_positions[0],
				   Y_positions[0],
				   X_positions[1],
				   Y_positions[1],
				   poly_lines[1]);

    introduce_PointInsideHalfPlane(z_solver,
				   X_positions[0],
				   Y_positions[0],
				   X_positions[1],
				   Y_positions[1],
				   poly_lines[2]);

    introduce_PointInsideHalfPlane(z_solver,
				   X_positions[0],
				   Y_positions[0],
				   X_positions[1],
				   Y_positions[1],
				   poly_lines[3]);        

    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    switch (z_solver.check())
    {
    case z3::sat:
    {
	printf("  SATISFIABLE\n");
	break;
    }
    case z3::unsat:	
    {
	printf("  UNSATISFIABLE\n");
	return;
	break;
    }
    case z3::unknown:
    {
	printf("  UNKNOWN\n");
	break;
    }
    default:
    {
	break;
    }
    }   

    z3::model z_model(z_solver.get_model());
    printf("Printing model:\n");
    cout << z_model << "\n";

    finish = clock();    

    printf("Printing interpretation:\n");    
    for (unsigned int i = 0; i < z_model.size(); ++i)
    {
	printf("Variable:%s  ", z_model[i].name().str().c_str());
	
	cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
	double value = z_model.get_const_interp(z_model[i]).as_double();
	
	printf("value: %.3f\n", value);
	
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 5 ... finished\n");    
}


Polygon polygon_1 = {{0, 0}, {50, 0}, {50, 50}, {0, 50}};
//Polygon polygon_1 = {{scale_(0), scale_(0)}, {scale_(50), scale_(0)}, {scale_(50), scale_(50)}, {scale_(0), scale_(50)}};


void test_polygon_6(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 6 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    
    for (int i = 0; i < poly_line_count; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < poly_line_count; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }
    
    z3::solver z_solver(z_context);

    introduce_PointOutsidePolygon(z_solver,
				  z_context,
				  X_positions[0],
				  Y_positions[0],
				  X_positions[1],
				  Y_positions[1],
				  polygon_1);

    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    switch (z_solver.check())
    {
    case z3::sat:
    {
	printf("  SATISFIABLE\n");
	break;
    }
    case z3::unsat:	
    {
	printf("  UNSATISFIABLE\n");
	return;
	break;
    }
    case z3::unknown:
    {
	printf("  UNKNOWN\n");
	break;
    }
    default:
    {
	break;
    }
    }   

    z3::model z_model(z_solver.get_model());
    printf("Printing model:\n");
    cout << z_model << "\n";

    finish = clock();    

    printf("Printing interpretation:\n");    
    for (unsigned int i = 0; i < z_model.size(); ++i)
    {
	printf("Variable:%s  ", z_model[i].name().str().c_str());
	
	cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
	double value = z_model.get_const_interp(z_model[i]).as_double();

	z3::expr valo_1 = z_model.get_const_interp(z_model[i]);
	z3::expr deco_1 = expr(z_context.real_const("deco_1"));

	z3::expr lino_1 = (valo_1 * deco_1 == 0);
	
	printf("value: %.3f\n", value);
	
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 6 ... finished\n");    
}


Polygon polygon_2 = {{0, 0}, {150, 0}, {150, 50}, {75, 120}, {0, 50} };
//Polygon polygon_2 = {{scale_(0), scale_(0)}, {scale_(150), scale_(0)}, {scale_(150), scale_(50)}, {scale_(75), scale_(120)}, {scale_(0), scale_(50)} };


void test_polygon_7(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 7 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T1_parameters(z_context);
    z3::expr_vector T2_parameters(z_context);    
    
    for (int i = 0; i < 2; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < 2; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (unsigned int i = 0; i < polygon_1.points.size(); ++i)
    {
	string name = "t1_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T1_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_2.points.size(); ++i)
    {
	string name = "t2_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T2_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }        
    
    z3::solver z_solver(z_context);

   introduce_DecisionBox(z_solver, X_positions[0], Y_positions[0], 200, 200);
   introduce_DecisionBox(z_solver, X_positions[1], Y_positions[1], 200, 200);    

    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[0],
				    Y_positions[0],
				    polygon_1,
				    X_positions[1],
				    Y_positions[1],				    
				    polygon_2);

    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    switch (z_solver.check())
    {
    case z3::sat:
    {
	printf("  SATISFIABLE\n");
	break;
    }
    case z3::unsat:	
    {
	printf("  UNSATISFIABLE\n");
	return;
	break;
    }
    case z3::unknown:
    {
	printf("  UNKNOWN\n");
	break;
    }
    default:
    {
	break;
    }
    }   

    z3::model z_model(z_solver.get_model());
    printf("Printing model:\n");
    cout << z_model << "\n";

    finish = clock();    

    double poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y;
    poly_1_pos_x = poly_1_pos_y = poly_2_pos_x = poly_2_pos_y = 0.0;    
    
    printf("Printing interpretation:\n");    
    for (unsigned int i = 0; i < z_model.size(); ++i)
    {
	printf("Variable:%s  ", z_model[i].name().str().c_str());
	
	cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
	double value = z_model.get_const_interp(z_model[i]).as_double();
	printf("value: %.3f\n", value);	

	if (z_model[i].name().str() == "x_pos-0")
	{
	    poly_1_pos_x = value;
	}
	else if (z_model[i].name().str() == "y_pos-0")
	{
	    poly_1_pos_y = value;
	}
	else if (z_model[i].name().str() == "x_pos-1")
	{
	    poly_2_pos_x = value;
	}
	else if (z_model[i].name().str() == "y_pos-1")
	{
	    poly_2_pos_y = value;
	}	
	
	
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }

    printf("Positions: %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y);

    /*
    for (int i = 0; i < 2; ++i)
    {	
	double value = X_positions[i].as_double();
	printf("Orig X: %.3f\n", value);

	value = Y_positions[i].as_double();
	printf("Orig Y: %.3f\n", value);	
    }
    */
    
    SVG preview_svg("polygon_test_7.svg");

//    preview_svg.draw(polygon_1);
//    preview_svg.draw(polygon_2);
    
    preview_svg.Close();    

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 7 ... finished\n");    
}


Polygon scale_UP(const Polygon &polygon)
{
    Polygon poly = polygon;

    for (unsigned int i = 0; i < poly.points.size(); ++i)
    {
	poly.points[i] = Point(poly.points[i].x() * SCALE_FACTOR, poly.points[i].y() * SCALE_FACTOR);
    }

    return poly;
}


Polygon scale_UP(const Polygon &polygon, double x_pos, double y_pos)
{
    Polygon poly = polygon;

    for (unsigned int i = 0; i < poly.points.size(); ++i)
    {
	poly.points[i] = Point(poly.points[i].x() * SCALE_FACTOR + x_pos * SCALE_FACTOR, poly.points[i].y() * SCALE_FACTOR + y_pos * SCALE_FACTOR);
    }

    return poly;    
}


Polygon polygon_3 = {{40, 0}, {80, 40}, {40, 80}, {0, 40}};
//Polygon polygon_3 = {{20, 0}, {40, 0}, {60, 30}, {30, 50}, {0, 30}};

void test_polygon_8(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 8 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T1_parameters(z_context);
    z3::expr_vector T2_parameters(z_context);
    z3::expr_vector T3_parameters(z_context);        
    
    for (int i = 0; i < 3; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < 3; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (unsigned int i = 0; i < polygon_1.points.size(); ++i)
    {
	string name = "t1_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T1_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_2.points.size(); ++i)
    {
	string name = "t2_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T2_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_3.points.size(); ++i)
    {
	string name = "t3_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T3_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }            
    
    z3::solver z_solver(z_context);

    /*
    introduce_DecisionBox(z_solver, X_positions[0], Y_positions[0], 200, 200);
    introduce_DecisionBox(z_solver, X_positions[1], Y_positions[1], 200, 200);    
    */

    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[0],
				    Y_positions[0],
				    polygon_1,
				    X_positions[1],
				    Y_positions[1],				    
				    polygon_2);

    /*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[1],
				    Y_positions[1],
				    polygon_2,
				    X_positions[0],
				    Y_positions[0],				    
				    polygon_1);    
    */

    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[0],
					 Y_positions[0],
					 polygon_1,
					 X_positions[1],
					 Y_positions[1],
					 polygon_2);
    
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[1],
				    Y_positions[1],
				    polygon_2,
				    X_positions[2],
				    Y_positions[2],				    
				    polygon_3);

    /*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[2],
				    Y_positions[2],
				    polygon_3,
				    X_positions[1],
				    Y_positions[1],				    
				    polygon_2);
    */

    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[1],
					 Y_positions[1],
					 polygon_2,
					 X_positions[2],
					 Y_positions[2],
					 polygon_3);

    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[0],
				    Y_positions[0],
				    polygon_1,
				    X_positions[2],
				    Y_positions[2],				    
				    polygon_3);

/*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[2],
				    Y_positions[2],
				    polygon_3,
				    X_positions[0],
				    Y_positions[0],				    
				    polygon_1);
*/

    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[0],
					 Y_positions[0],
					 polygon_1,
					 X_positions[2],
					 Y_positions[2],
					 polygon_3);	

  
    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    int last_solvable_decision_box_size = -1;
    double poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y;
    poly_1_pos_x = poly_1_pos_y = poly_2_pos_x = poly_2_pos_y = poly_3_pos_x = poly_3_pos_y = 0.0;    
    
    for (int decision_box_size = 300; decision_box_size > 10; decision_box_size -= 4)
    {
	z3::expr_vector decision_box_assumptions(z_context);
	
	assume_DecisionBox(X_positions[0], Y_positions[0], decision_box_size, decision_box_size, decision_box_assumptions);
	assume_DecisionBox(X_positions[1], Y_positions[1], decision_box_size, decision_box_size, decision_box_assumptions);
	assume_DecisionBox(X_positions[2], Y_positions[2], decision_box_size, decision_box_size, decision_box_assumptions);	

	bool sat = false;
	
	switch (z_solver.check(decision_box_assumptions))
	{
	case z3::sat:
	{	    
	    printf("  SATISFIABLE\n");
	    last_solvable_decision_box_size = decision_box_size;
	    sat = true;	    
	    break;
	}
	case z3::unsat:	
	{
	    printf("  UNSATISFIABLE\n");
	    sat = false;	    
	    break;
	}
	case z3::unknown:
	{
	    printf("  UNKNOWN\n");
	    break;
	}
	default:
	{
	    break;
	}
	}

	if (sat)
	{
	    z3::model z_model(z_solver.get_model());
	    printf("Printing model:\n");
	    cout << z_model << "\n";
    
	    printf("Printing interpretation:\n");    
	    for (unsigned int i = 0; i < z_model.size(); ++i)
	    {
		printf("Variable:%s  ", z_model[i].name().str().c_str());
		
		cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
		double value = z_model.get_const_interp(z_model[i]).as_double();
		printf("value: %.3f\n", value);	
	    
		if (z_model[i].name().str() == "x_pos-0")
		{
		    poly_1_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-0")
		{
		    poly_1_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-1")
		{
		    poly_2_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-1")
		{
		    poly_2_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-2")
		{
		    poly_3_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-2")
		{
		    poly_3_pos_y = value;
		}		
	    }
	}
	else
	{
	    break;
	}
		
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }
    finish = clock();

    printf("Solvable decision box: %d\n", last_solvable_decision_box_size);
    printf("Positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y);

    /*
    for (int i = 0; i < 2; ++i)
    {	
	double value = X_positions[i].as_double();
	printf("Orig X: %.3f\n", value);

	value = Y_positions[i].as_double();
	printf("Orig Y: %.3f\n", value);	
    }
    */
    
    SVG preview_svg("polygon_test_8.svg");

    Polygon display_polygon_1 = scale_UP(polygon_1, poly_1_pos_x, poly_1_pos_y);
    Polygon display_polygon_2 = scale_UP(polygon_2, poly_2_pos_x, poly_2_pos_y);
    Polygon display_polygon_3 = scale_UP(polygon_3, poly_3_pos_x, poly_3_pos_y);    

    preview_svg.draw(display_polygon_1, "green");
    preview_svg.draw(display_polygon_2, "blue");
    preview_svg.draw(display_polygon_3, "red");    
    
    preview_svg.Close();    

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 8 ... finished\n");    
}


void test_polygon_9(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 9 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T1_parameters(z_context);
    z3::expr_vector T2_parameters(z_context);
    z3::expr_vector T3_parameters(z_context);        
    
    for (int i = 0; i < 3; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < 3; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (unsigned int i = 0; i < polygon_1.points.size(); ++i)
    {
	string name = "t1_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T1_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_2.points.size(); ++i)
    {
	string name = "t2_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T2_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_3.points.size(); ++i)
    {
	string name = "t3_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T3_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }            
    
    z3::solver z_solver(z_context);

    /*
    introduce_DecisionBox(z_solver, X_positions[0], Y_positions[0], 200, 200);
    introduce_DecisionBox(z_solver, X_positions[1], Y_positions[1], 200, 200);    
    */

    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[0],
				    Y_positions[0],
				    polygon_1,
				    X_positions[1],
				    Y_positions[1],				    
				    polygon_2);
/*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[1],
				    Y_positions[1],
				    polygon_2,
				    X_positions[0],
				    Y_positions[0],				    
				    polygon_1);    
*/

    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[0],
					 Y_positions[0],
					 polygon_1,
					 X_positions[1],
					 Y_positions[1],
					 polygon_2);
    
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[1],
				    Y_positions[1],
				    polygon_2,
				    X_positions[2],
				    Y_positions[2],				    
				    polygon_3);
/*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[2],
				    Y_positions[2],
				    polygon_3,
				    X_positions[1],
				    Y_positions[1],				    
				    polygon_2);
*/
    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[1],
					 Y_positions[1],
					 polygon_2,
					 X_positions[2],
					 Y_positions[2],
					 polygon_3);

    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[0],
				    Y_positions[0],
				    polygon_1,
				    X_positions[2],
				    Y_positions[2],				    
				    polygon_3);
/*
    introduce_PolygonOutsidePolygon(z_solver,
				    z_context,
				    X_positions[2],
				    Y_positions[2],
				    polygon_3,
				    X_positions[0],
				    Y_positions[0],				    
				    polygon_1);
*/

    introduce_PolygonLineNonIntersection(z_solver,
					 z_context,
					 X_positions[0],
					 Y_positions[0],
					 polygon_1,
					 X_positions[2],
					 Y_positions[2],
					 polygon_3);	

  
    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    int last_solvable_bounding_box_size = -1;
    double poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y;
    poly_1_pos_x = poly_1_pos_y = poly_2_pos_x = poly_2_pos_y = poly_3_pos_x = poly_3_pos_y = 0.0;
    
    for (int bounding_box_size = 300; bounding_box_size > 10; bounding_box_size -= 4)
    {
	z3::expr_vector bounding_box_assumptions(z_context);
	
	assume_BedBoundingBox(X_positions[0], Y_positions[0], polygon_1, bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[1], Y_positions[1], polygon_2, bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[2], Y_positions[2], polygon_3, bounding_box_size, bounding_box_size, bounding_box_assumptions);	

	bool sat = false;
	
	switch (z_solver.check(bounding_box_assumptions))
	{
	case z3::sat:
	{	    
	    printf("  SATISFIABLE\n");
	    last_solvable_bounding_box_size = bounding_box_size;
	    sat = true;	    
	    break;
	}
	case z3::unsat:	
	{
	    printf("  UNSATISFIABLE\n");
	    sat = false;	    
	    break;
	}
	case z3::unknown:
	{
	    printf("  UNKNOWN\n");
	    break;
	}
	default:
	{
	    break;
	}
	}

	if (sat)
	{
	    z3::model z_model(z_solver.get_model());
	    printf("Printing model:\n");
	    cout << z_model << "\n";
    
	    printf("Printing interpretation:\n");    
	    for (unsigned int i = 0; i < z_model.size(); ++i)
	    {
		printf("Variable:%s  ", z_model[i].name().str().c_str());
		
		cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
		double value = z_model.get_const_interp(z_model[i]).as_double();
		printf("value: %.3f\n", value);	
	    
		if (z_model[i].name().str() == "x_pos-0")
		{
		    poly_1_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-0")
		{
		    poly_1_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-1")
		{
		    poly_2_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-1")
		{
		    poly_2_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-2")
		{
		    poly_3_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-2")
		{
		    poly_3_pos_y = value;
		}		
	    }
	}
	else
	{
	    break;
	}
		
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }
    finish = clock();

    printf("Solvable bounding box: %d\n", last_solvable_bounding_box_size);
    printf("Positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y);

    /*
    for (int i = 0; i < 2; ++i)
    {	
	double value = X_positions[i].as_double();
	printf("Orig X: %.3f\n", value);

	value = Y_positions[i].as_double();
	printf("Orig Y: %.3f\n", value);	
    }
    */
    
    SVG preview_svg("polygon_test_9.svg");

    Polygon display_polygon_1 = scale_UP(polygon_1, poly_1_pos_x, poly_1_pos_y);
    Polygon display_polygon_2 = scale_UP(polygon_2, poly_2_pos_x, poly_2_pos_y);
    Polygon display_polygon_3 = scale_UP(polygon_3, poly_3_pos_x, poly_3_pos_y);    

    preview_svg.draw(display_polygon_1, "green");
    preview_svg.draw(display_polygon_2, "blue");
    preview_svg.draw(display_polygon_3, "red");    
    
    preview_svg.Close();    

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 9 ... finished\n");    
}


Polygon polygon_4 = {{20, 0}, {40, 0}, {60, 30}, {30, 50}, {0, 30}};

void test_polygon_10(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 10 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T1_parameters(z_context);
    z3::expr_vector T2_parameters(z_context);
    z3::expr_vector T3_parameters(z_context);
    z3::expr_vector T4_parameters(z_context);            
    
    for (int i = 0; i < 4; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < 4; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (unsigned int i = 0; i < polygon_1.points.size(); ++i)
    {
	string name = "t1_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T1_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_2.points.size(); ++i)
    {
	string name = "t2_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T2_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_3.points.size(); ++i)
    {
	string name = "t3_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T3_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_4.points.size(); ++i)
    {
	string name = "t4_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T4_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }                
    
    z3::solver z_solver(z_context);

    vector<Polygon> polygons;
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);    
    
    /*
    introduce_DecisionBox(z_solver, X_positions[0], Y_positions[0], 200, 200);
    introduce_DecisionBox(z_solver, X_positions[1], Y_positions[1], 200, 200);    
    */

    /*
    for (int i = 0; i < 3; ++i)
    {
	for (int j = i + 1; j < 4; ++j)
	{	    
	    introduce_PolygonOutsidePolygon(z_solver,
					    z_context,
					    X_positions[i],
					    Y_positions[i],
					    polygons[i],
					    X_positions[j],
					    Y_positions[j],				    
					    polygons[j]);
	}
    }
    */

/*
        introduce_PolygonOutsidePolygon(z_solver,
					z_context,
					X_positions[0],
					Y_positions[0],
					polygons[0],
					X_positions[2],
					Y_positions[2],				    
					polygons[2]);    

	
	introduce_PolygonOutsidePolygon(z_solver,
					z_context,
					X_positions[0],
					Y_positions[0],
					polygons[0],
					X_positions[1],
					Y_positions[1],				    
					polygons[1]);
	
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[0],
				       Y_positions[0],
				       polygons[0],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);
       
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[1],
				       Y_positions[1],
				       polygons[1],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);
       
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[1],
				       Y_positions[1],
				       polygons[1],
				       X_positions[2],
				       Y_positions[2],				    
				       polygons[2]);

       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[2],
				       Y_positions[2],
				       polygons[2],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);       
*/
       
/*
    introduce_PolygonWeakNonoverlapping(z_solver,
					z_context,
					X_positions,
					Y_positions,
					polygons);
*/   
    introduce_PolygonStrongNonoverlapping(z_solver,
					  z_context,
					  X_positions,
					  Y_positions,
					  polygons);    
  
    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    int last_solvable_bounding_box_size = -1;
    double poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y, poly_4_pos_x, poly_4_pos_y;
    poly_1_pos_x = poly_1_pos_y = poly_2_pos_x = poly_2_pos_y = poly_3_pos_x = poly_3_pos_y = poly_4_pos_x = poly_4_pos_y = 0.0;    
   
    for (int bounding_box_size = 300; bounding_box_size > 10; bounding_box_size -= 4)
    {
	z3::expr_vector bounding_box_assumptions(z_context);

	//assume_BedBoundingBox(X_positions, Y_positions, polygons, bounding_box_size, bounding_box_size, bounding_box_assumptions);	
	assume_BedBoundingBox(X_positions[0], Y_positions[0], polygons[0], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[1], Y_positions[1], polygons[1], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[2], Y_positions[2], polygons[2], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[3], Y_positions[3], polygons[3], bounding_box_size, bounding_box_size, bounding_box_assumptions);		
	
	bool sat = false;
	
	switch (z_solver.check(bounding_box_assumptions))
	{
	case z3::sat:
	{	    
	    printf("  SATISFIABLE\n");
	    last_solvable_bounding_box_size = bounding_box_size;
	    sat = true;	    
	    break;
	}
	case z3::unsat:	
	{
	    printf("  UNSATISFIABLE\n");
	    sat = false;	    
	    break;
	}
	case z3::unknown:
	{
	    printf("  UNKNOWN\n");
	    break;
	}
	default:
	{
	    break;
	}
	}

	if (sat)
	{
	    z3::model z_model(z_solver.get_model());
	    printf("Printing model:\n");
	    cout << z_model << "\n";
    
	    printf("Printing interpretation:\n");    
	    for (unsigned int i = 0; i < z_model.size(); ++i)
	    {
		printf("Variable:%s  ", z_model[i].name().str().c_str());
		
		cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
		double value = z_model.get_const_interp(z_model[i]).as_double();
		printf("value: %.3f\n", value);	
	    
		if (z_model[i].name().str() == "x_pos-0")
		{
		    poly_1_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-0")
		{
		    poly_1_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-1")
		{
		    poly_2_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-1")
		{
		    poly_2_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-2")
		{
		    poly_3_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-2")
		{
		    poly_3_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-3")
		{
		    poly_4_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-3")
		{
		    poly_4_pos_y = value;
		}				
	    }
	}
	else
	{
	    break;
	}	
		
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }
    finish = clock();

    printf("Solvable bounding box: %d\n", last_solvable_bounding_box_size);
    printf("Positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x,
	                                                                  poly_1_pos_y,
                                                                  	  poly_2_pos_x,
                                                                  	  poly_2_pos_y,
                                                                  	  poly_3_pos_x,
                                                                  	  poly_3_pos_y,
                                                                  	  poly_4_pos_x,
                                                                  	  poly_4_pos_y);

    /*
    for (int i = 0; i < 2; ++i)
    {	
	double value = X_positions[i].as_double();
	printf("Orig X: %.3f\n", value);

	value = Y_positions[i].as_double();
	printf("Orig Y: %.3f\n", value);	
    }
    */
    
    SVG preview_svg("polygon_test_10.svg");

    Polygon display_polygon_1 = scale_UP(polygons[0], poly_1_pos_x, poly_1_pos_y);
    Polygon display_polygon_2 = scale_UP(polygons[1], poly_2_pos_x, poly_2_pos_y);
    Polygon display_polygon_3 = scale_UP(polygons[2], poly_3_pos_x, poly_3_pos_y);
    Polygon display_polygon_4 = scale_UP(polygons[3], poly_4_pos_x, poly_4_pos_y);        

    preview_svg.draw(display_polygon_1, "green");
    preview_svg.draw(display_polygon_2, "blue");
    preview_svg.draw(display_polygon_3, "red");
    preview_svg.draw(display_polygon_4, "grey");        
    
    preview_svg.Close();    

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 10 ... finished\n");    
}


void test_polygon_11(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 11 ...\n");

    start = clock();

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);
    z3::expr_vector T1_parameters(z_context);
    z3::expr_vector T2_parameters(z_context);
    z3::expr_vector T3_parameters(z_context);
    z3::expr_vector T4_parameters(z_context);            
    
    for (int i = 0; i < 4; ++i)
    {
	printf("i:%d\n", i);
	string name = "x_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	X_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (int i = 0; i < 4; ++i)
    {
	string name = "y_pos-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	Y_positions.push_back(expr(z_context.real_const(name.c_str())));
    }

    for (unsigned int i = 0; i < polygon_1.points.size(); ++i)
    {
	string name = "t1_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T1_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_2.points.size(); ++i)
    {
	string name = "t2_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T2_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_3.points.size(); ++i)
    {
	string name = "t3_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T3_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }

    for (unsigned int i = 0; i < polygon_4.points.size(); ++i)
    {
	string name = "t4_par-" + to_string(i);
	printf("name: %s\n", name.c_str());
	
	T4_parameters.push_back(expr(z_context.real_const(name.c_str())));	
    }                
    
    z3::solver z_solver(z_context);

    vector<Polygon> polygons;
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);    
    
    /*
    introduce_DecisionBox(z_solver, X_positions[0], Y_positions[0], 200, 200);
    introduce_DecisionBox(z_solver, X_positions[1], Y_positions[1], 200, 200);    
    */

    /*
    for (int i = 0; i < 3; ++i)
    {
	for (int j = i + 1; j < 4; ++j)
	{	    
	    introduce_PolygonOutsidePolygon(z_solver,
					    z_context,
					    X_positions[i],
					    Y_positions[i],
					    polygons[i],
					    X_positions[j],
					    Y_positions[j],				    
					    polygons[j]);
	}
    }
    */

/*
        introduce_PolygonOutsidePolygon(z_solver,
					z_context,
					X_positions[0],
					Y_positions[0],
					polygons[0],
					X_positions[2],
					Y_positions[2],				    
					polygons[2]);    

	
	introduce_PolygonOutsidePolygon(z_solver,
					z_context,
					X_positions[0],
					Y_positions[0],
					polygons[0],
					X_positions[1],
					Y_positions[1],				    
					polygons[1]);
	
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[0],
				       Y_positions[0],
				       polygons[0],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);
       
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[1],
				       Y_positions[1],
				       polygons[1],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);
       
       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[1],
				       Y_positions[1],
				       polygons[1],
				       X_positions[2],
				       Y_positions[2],				    
				       polygons[2]);

       introduce_PolygonOutsidePolygon(z_solver,
				       z_context,
				       X_positions[2],
				       Y_positions[2],
				       polygons[2],
				       X_positions[3],
				       Y_positions[3],				    
				       polygons[3]);       
*/
       
    introduce_PolygonWeakNonoverlapping(z_solver,
					z_context,
					X_positions,
					Y_positions,
					polygons);
  
    printf("Printing solver status:\n");
    cout << z_solver << "\n";
    
    printf("Printing smt status:\n");
    cout << z_solver.to_smt2() << "\n";

    int last_solvable_bounding_box_size = -1;
    double poly_1_pos_x, poly_1_pos_y, poly_2_pos_x, poly_2_pos_y, poly_3_pos_x, poly_3_pos_y, poly_4_pos_x, poly_4_pos_y;
    poly_1_pos_x = poly_1_pos_y = poly_2_pos_x = poly_2_pos_y = poly_3_pos_x = poly_3_pos_y = poly_4_pos_x = poly_4_pos_y = 0.0;    
   
    for (int bounding_box_size = 200; bounding_box_size > 10; bounding_box_size -= 4)
    {
	printf("BB: %d\n", bounding_box_size);
	z3::expr_vector bounding_box_assumptions(z_context);

	//assume_BedBoundingBox(X_positions, Y_positions, polygons, bounding_box_size, bounding_box_size, bounding_box_assumptions);	
	assume_BedBoundingBox(X_positions[0], Y_positions[0], polygons[0], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[1], Y_positions[1], polygons[1], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[2], Y_positions[2], polygons[2], bounding_box_size, bounding_box_size, bounding_box_assumptions);
	assume_BedBoundingBox(X_positions[3], Y_positions[3], polygons[3], bounding_box_size, bounding_box_size, bounding_box_assumptions);		
	
	bool sat = false;
	
	switch (z_solver.check(bounding_box_assumptions))
	{
	case z3::sat:
	{	    
	    printf("  SATISFIABLE\n");
	    sat = true;	    
	    break;
	}
	case z3::unsat:	
	{
	    printf("  UNSATISFIABLE\n");
	    sat = false;	    
	    break;
	}
	case z3::unknown:
	{
	    printf("  UNKNOWN\n");
	    break;
	}
	default:
	{
	    break;
	}
	}

	if (sat)
	{
	    z3::model z_model(z_solver.get_model());
	    printf("Printing model:\n");
	    cout << z_model << "\n";
    
	    printf("Printing interpretation:\n");    
	    for (unsigned int i = 0; i < z_model.size(); ++i)
	    {
		printf("Variable:%s  ", z_model[i].name().str().c_str());
		
		cout << z_model.get_const_interp(z_model[i]).as_double() << "\n";
		double value = z_model.get_const_interp(z_model[i]).as_double();
		printf("value: %.3f\n", value);	
	    
		if (z_model[i].name().str() == "x_pos-0")
		{
		    poly_1_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-0")
		{
		    poly_1_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-1")
		{
		    poly_2_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-1")
		{
		    poly_2_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-2")
		{
		    poly_3_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-2")
		{
		    poly_3_pos_y = value;
		}
		else if (z_model[i].name().str() == "x_pos-3")
		{
		    poly_4_pos_x = value;
		}
		else if (z_model[i].name().str() == "y_pos-3")
		{
		    poly_4_pos_y = value;
		}				
	    }

	    printf("preRefined positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x,
	                                                                  poly_1_pos_y,
                                                                  	  poly_2_pos_x,
                                                                  	  poly_2_pos_y,
                                                                  	  poly_3_pos_x,
                                                                  	  poly_3_pos_y,
                                                                  	  poly_4_pos_x,
                                                                  	  poly_4_pos_y);				    

	    while (true)
	    {
		vector<double> dec_values_X;
		dec_values_X.push_back(poly_1_pos_x);
		dec_values_X.push_back(poly_2_pos_x);
		dec_values_X.push_back(poly_3_pos_x);
		dec_values_X.push_back(poly_4_pos_x);
	    
		vector<double> dec_values_Y;
		dec_values_Y.push_back(poly_1_pos_y);
		dec_values_Y.push_back(poly_2_pos_y);
		dec_values_Y.push_back(poly_3_pos_y);
		dec_values_Y.push_back(poly_4_pos_y);	    
		
		bool refined = refine_PolygonWeakNonoverlapping(z_solver,
								z_context,
								X_positions,
								Y_positions,
								dec_values_X,
								dec_values_Y,
								polygons);

		bool refined_sat = false;

		if (refined)
		{
		    switch (z_solver.check(bounding_box_assumptions))
		    {
		    case z3::sat:
		    {	    
			printf("  sat\n");
			refined_sat = true;	    
			break;
		    }
		    case z3::unsat:	
		    {
			printf("  unsat\n");
			refined_sat = false;	    
			break;
		    }
		    case z3::unknown:
		    {
			printf("  unknown\n");
			break;
		    }
		    default:
		    {
			break;
		    }
		    }

		    if (refined_sat)
		    {
			z3::model z_model(z_solver.get_model());			
			printf("Printing model:\n");
			cout << z_model << "\n";
    
			for (unsigned int i = 0; i < z_model.size(); ++i)
			{
			    //printf("Variable:%s  ", z_model[i].name().str().c_str());						    
			    double value = z_model.get_const_interp(z_model[i]).as_double();
			
			    if (z_model[i].name().str() == "x_pos-0")
			    {
				poly_1_pos_x = value;
			    }
			    else if (z_model[i].name().str() == "y_pos-0")
			    {
				poly_1_pos_y = value;
			    }
			    else if (z_model[i].name().str() == "x_pos-1")
			    {
				poly_2_pos_x = value;
			    }
			    else if (z_model[i].name().str() == "y_pos-1")
			    {
				poly_2_pos_y = value;
			    }
			    else if (z_model[i].name().str() == "x_pos-2")
			    {
				poly_3_pos_x = value;
			    }
			    else if (z_model[i].name().str() == "y_pos-2")
			    {
				poly_3_pos_y = value;
			    }
			    else if (z_model[i].name().str() == "x_pos-3")
			    {
				poly_4_pos_x = value;
			    }
			    else if (z_model[i].name().str() == "y_pos-3")
			    {
				poly_4_pos_y = value;
			    }
			}
			printf("Refined positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x,
	                                                                  poly_1_pos_y,
                                                                  	  poly_2_pos_x,
                                                                  	  poly_2_pos_y,
                                                                  	  poly_3_pos_x,
                                                                  	  poly_3_pos_y,
                                                                  	  poly_4_pos_x,
                                                                  	  poly_4_pos_y);			
		    }
		    else
		    {
			break;
		    }
		}
		else
		{
		    last_solvable_bounding_box_size = bounding_box_size;
		    break;
		}
	    }
	}
	else
	{
	    break;
	}	
		
	//cout << float(z_model[i]) << "\n";
        /*
	switch (z_model.get_const_interp(z_model[i]).bool_value())
	{
	case Z3_L_FALSE:
	{
	    printf("   value: FALSE\n");
	    break;
	}
	case Z3_L_TRUE:
	{
	    printf("   value: TRUE\n");
	    break;
	}
	case Z3_L_UNDEF:
	{
	    printf("   value: UNDEF\n");
		break;
	}	    
	default:
	{
		break;
	}
	}
	*/
    }
    finish = clock();

    printf("Solvable bounding box: %d\n", last_solvable_bounding_box_size);
    printf("Positions: %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n", poly_1_pos_x,
	                                                                  poly_1_pos_y,
                                                                  	  poly_2_pos_x,
                                                                  	  poly_2_pos_y,
                                                                  	  poly_3_pos_x,
                                                                  	  poly_3_pos_y,
                                                                  	  poly_4_pos_x,
                                                                  	  poly_4_pos_y);

    /*
    for (int i = 0; i < 2; ++i)
    {	
	double value = X_positions[i].as_double();
	printf("Orig X: %.3f\n", value);

	value = Y_positions[i].as_double();
	printf("Orig Y: %.3f\n", value);	
    }
    */
    
    SVG preview_svg("polygon_test_11.svg");

    Polygon display_polygon_1 = scale_UP(polygons[0], poly_1_pos_x, poly_1_pos_y);
    Polygon display_polygon_2 = scale_UP(polygons[1], poly_2_pos_x, poly_2_pos_y);
    Polygon display_polygon_3 = scale_UP(polygons[2], poly_3_pos_x, poly_3_pos_y);
    Polygon display_polygon_4 = scale_UP(polygons[3], poly_4_pos_x, poly_4_pos_y);        

    preview_svg.draw(display_polygon_1, "green");
    preview_svg.draw(display_polygon_2, "blue");
    preview_svg.draw(display_polygon_3, "red");
    preview_svg.draw(display_polygon_4, "grey");        
    
    preview_svg.Close();    

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 11 ... finished\n");    
}


void test_polygon_12(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 12 ...\n");

    start = clock();

    SolverConfiguration solver_configuration;

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);

    std::vector<double> X_values;
    std::vector<double> Y_values;
    
    string_map dec_var_names_map;
       
    z3::solver z_solver(z_context);

    vector<Polygon> polygons;
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);    
    
    build_WeakPolygonNonoverlapping(z_solver, z_context, polygons, X_positions, Y_positions, X_values, Y_values, dec_var_names_map);

    bool optimized = optimize_WeakPolygonNonoverlapping(z_solver,
							z_context,
							solver_configuration,
							X_positions,
							Y_positions,
							X_values,
							Y_values,							
							dec_var_names_map,
							polygons);

    finish = clock();

    if (optimized)
    {
	printf("Polygon positions:\n");
	for (unsigned int i = 0; i < polygons.size(); ++i)
	{
	    printf("  %.3f, %.3f\n", X_values[i], Y_values[i]);
	}
    
	SVG preview_svg("polygon_test_12.svg");
	
	for (unsigned int i = 0; i < polygons.size(); ++i)
	{
	    Polygon display_polygon = scale_UP(polygons[i], X_values[i], Y_values[i]);
	    
	    string color;
	    
	    switch(i)
	    {
	    case 0:
	    {
		color = "green";
		break;
	    }
	    case 1:
	    {
		color = "blue";
		break;
	    }
	    case 2:
	    {
		color = "red";	    
		break;
	    }
	    case 3:
	    {
		color = "grey";	    
		break;
	    }
	    case 4:
	    {
		color = "cyan";
		break;
	    }
	    case 5:
	    {
		color = "magenta";
		break;
	    }		
	    default:
	    {
		break;
	    }
	    }

	    preview_svg.draw(display_polygon, color);
	}
    
	preview_svg.Close();
    }
    else
    {
	printf("Polygon optimization FAILED.\n");
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 12 ... finished\n");    
}


void test_polygon_13(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 13 ...\n");

    start = clock();

    SolverConfiguration solver_configuration;

    z3::context z_context;    
    z3::expr_vector X_positions(z_context);
    z3::expr_vector Y_positions(z_context);

    std::vector<double> X_values;
    std::vector<double> Y_values;
    
    string_map dec_var_names_map;

    /*
    z3::params z_parameters(z_context);
    z_parameters.set("timeout", 1000);
    */
    Z3_global_param_set("timeout", "8000");
       
    z3::solver z_solver(z_context);

    vector<Polygon> polygons;
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);
    
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);    

    build_WeakPolygonNonoverlapping(z_solver, z_context, polygons, X_positions, Y_positions, X_values, Y_values, dec_var_names_map);

    bool optimized = optimize_WeakPolygonNonoverlapping(z_solver,
							z_context,
							solver_configuration,
							X_positions,
							Y_positions,
							X_values,
							Y_values,							
							dec_var_names_map,
							polygons);

    finish = clock();

    if (optimized)
    {
	printf("Polygon positions:\n");
	for (unsigned int i = 0; i < polygons.size(); ++i)
	{
	    printf("  %.3f, %.3f\n", X_values[i], Y_values[i]);
	}
    
	SVG preview_svg("polygon_test_13.svg");
	
	for (unsigned int i = 0; i < polygons.size(); ++i)
	{
	    Polygon display_polygon = scale_UP(polygons[i], X_values[i], Y_values[i]);
	    
	    string color;
	    
	    switch(i)
	    {
	    case 0:
	    {
		color = "green";
		break;
	    }
	    case 1:
	    {
		color = "blue";
		break;
	    }
	    case 2:
	    {
		color = "red";	    
		break;
	    }
	    case 3:
	    {
		color = "grey";	    
		break;
	    }
	    case 4:
	    {
		color = "cyan";
		break;
	    }
	    case 5:
	    {
		color = "magenta";
		break;
	    }
	    case 6:
	    {
		color = "yellow";
		break;
	    }
	    case 7:
	    {
		color = "black";
		break;
	    }
	    case 8:
	    {
		color = "indigo";
		break;
	    }
	    case 9:
	    {
		color = "olive";
		break;
	    }
	    case 10:
	    {
		color = "aqua";
		break;
	    }
	    case 11:
	    {
		color = "violet";
		break;
	    }			    	    	    
	    default:
	    {
		break;
	    }
	    }

	    preview_svg.draw(display_polygon, color);
	}
    
	preview_svg.Close();
    }
    else
    {
	printf("Polygon optimization FAILED.\n");
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 13 ... finished\n");    
}


void test_polygon_14(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 14 ...\n");

    start = clock();

    SolverConfiguration solver_configuration;

    vector<Polygon> polygons;
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);
    
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);    
    
    vector<int> decided;
    vector<int> undecided;

    vector<Rational> poly_positions_X;
    vector<Rational> poly_positions_Y;
    poly_positions_X.resize(polygons.size());
    poly_positions_Y.resize(polygons.size());    

    bool optimized;    
    {
	z3::context z_context;    
	z3::expr_vector X_positions(z_context);
	z3::expr_vector Y_positions(z_context);
	
	vector<Rational> X_values;
	vector<Rational> Y_values;
	
	string_map dec_var_names_map;
    
    /*
    Z3_global_param_set("timeout", "8000");
    */
       
	z3::solver z_solver(z_context);
	
	X_values.resize(polygons.size());
	Y_values.resize(polygons.size());
	
	undecided.push_back(0);
	undecided.push_back(1);
	undecided.push_back(2);
	undecided.push_back(3);    
	
	build_WeakPolygonNonoverlapping(z_solver,
					z_context,
					polygons,
					X_positions,
					Y_positions,
					X_values,
					Y_values,
					decided,
					undecided,				    
					dec_var_names_map);
	
	optimized = optimize_WeakPolygonNonoverlapping(z_solver,
						       z_context,
						       solver_configuration,
						       X_positions,
						       Y_positions,
						       X_values,
						       Y_values,
						       decided,
						       undecided,						   
						       dec_var_names_map,
						       polygons);

	for (unsigned int i = 0; i < undecided.size(); ++i)
	{
	    poly_positions_X[undecided[i]] = X_values[undecided[i]];
	    poly_positions_Y[undecided[i]] = Y_values[undecided[i]];
	}

	printf("Optimized 1: %d\n", optimized);
    }
    
    {
	z3::context z_context;    
	z3::expr_vector X_positions(z_context);
	z3::expr_vector Y_positions(z_context);
	
	vector<Rational> X_values;
	vector<Rational> Y_values;
	
	string_map dec_var_names_map;
           
	z3::solver z_solver(z_context);
	
	X_values.resize(polygons.size());
	Y_values.resize(polygons.size());

	decided.push_back(0);
	decided.push_back(1);
	decided.push_back(2);
	decided.push_back(3);	
	
	for (unsigned int i = 0; i < decided.size(); ++i)
	{
	    X_values[decided[i]] = poly_positions_X[decided[i]];
	    Y_values[decided[i]] = poly_positions_Y[decided[i]];	    
	}
		
	undecided.clear();
	undecided.push_back(4);
	undecided.push_back(5);
	undecided.push_back(6);
	undecided.push_back(7);

	build_WeakPolygonNonoverlapping(z_solver,
					z_context,
					polygons,
					X_positions,
					Y_positions,
					X_values,
					Y_values,
					decided,
					undecided,				    
					dec_var_names_map);    
	
	optimized = optimize_WeakPolygonNonoverlapping(z_solver,
						       z_context,
						       solver_configuration,
						       X_positions,
						       Y_positions,
						       X_values,
						       Y_values,
						       decided,
						       undecided,
						       dec_var_names_map,
						       polygons);
	printf("Optimized 2: %d\n", optimized);
	
	decided.push_back(4);
	decided.push_back(5);
	decided.push_back(6);
	decided.push_back(7);       
	
	finish = clock();
	
	if (optimized)
	{
	    printf("Polygon positions:\n");
	    for (unsigned int i = 0; i < decided.size(); ++i)
	    {
		printf("  %.3f, %.3f\n", X_values[decided[i]].as_double(), Y_values[decided[i]].as_double());
	    }
	    
	    SVG preview_svg("polygon_test_14.svg");
	
	    for (unsigned int i = 0; i < decided.size(); ++i)
	    {
		Polygon display_polygon = scale_UP(polygons[decided[i]], X_values[decided[i]].as_double(), Y_values[decided[i]].as_double());
		
		string color;
		
		switch(i)
		{
		case 0:
		{
		    color = "green";
		    break;
		}
		case 1:
		{
		    color = "blue";
		    break;
		}
		case 2:
		{
		    color = "red";	    
		    break;
		}
		case 3:
		{
		    color = "grey";	    
		    break;
		}
		case 4:
		{
		    color = "cyan";
		    break;
		}
		case 5:
		{
		    color = "magenta";
		    break;
		}
		case 6:
		{
		    color = "yellow";
		    break;
		}
		case 7:
		{
		    color = "black";
		    break;
		}
		case 8:
		{
		    color = "indigo";
		    break;
		}
		case 9:
		{
		    color = "olive";
		    break;
		}
		case 10:
		{
		    color = "aqua";
		    break;
		}
		case 11:
		{
		    color = "violet";
		    break;
		}			    	    	    
		default:
		{
		    break;
		}
		}

		preview_svg.draw(display_polygon, color);
	    }
    
	    preview_svg.Close();
	}
	else
	{
	    printf("Polygon optimization FAILED.\n");
	}
    }

    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 14 ... finished\n");    
}


void test_polygon_15(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 15 ...\n");

    start = clock();

    SolverConfiguration solver_configuration;

    vector<Polygon> polygons;
    vector<int> remaining_polygons;
    vector<int> polygon_index_map;
    vector<int> decided_polygons;
    
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);

    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);
    
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    
    /*
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);                
    */
    for (unsigned int index = 0; index < polygons.size(); ++index)
    {
	polygon_index_map.push_back(index);
    }
    
    vector<Rational> poly_positions_X;
    vector<Rational> poly_positions_Y;
    
    /*
    poly_positions_X.resize(polygons.size());
    poly_positions_Y.resize(polygons.size());    
    */

    do
    {
	decided_polygons.clear();
	remaining_polygons.clear();
	
	bool optimized = optimize_SubglobalPolygonNonoverlapping(solver_configuration,
								 poly_positions_X,
								 poly_positions_Y,
								 polygons,
								 polygon_index_map,
								 decided_polygons,
								 remaining_polygons);
	
	if (optimized)
	{
	    printf("Polygon positions:\n");
	    for (unsigned int i = 0; i < decided_polygons.size(); ++i)
	    {
		printf("  %.3f, %.3f\n", poly_positions_X[decided_polygons[i]].as_double(), poly_positions_Y[decided_polygons[i]].as_double());
	    }
	    printf("Remaining polygons: %ld\n", remaining_polygons.size());
	    for (unsigned int i = 0; i < remaining_polygons.size(); ++i)
	    {
		printf("  %d\n", remaining_polygons[i]);
	    }
	
	    SVG preview_svg("polygon_test_15.svg");
	
	    for (unsigned int i = 0; i < decided_polygons.size(); ++i)
	    {
		Polygon display_polygon = scale_UP(polygons[decided_polygons[i]],
						   poly_positions_X[decided_polygons[i]].as_double(),
						   poly_positions_Y[decided_polygons[i]].as_double());
		
		string color;
		
		switch(i)
		{
		case 0:
		{
		    color = "green";
		    break;
		}
		case 1:
		{
		    color = "blue";
		    break;
		}
		case 2:
		{
		    color = "red";	    
		    break;
		}
		case 3:
		{
		    color = "grey";	    
		    break;
		}
		case 4:
		{
		    color = "cyan";
		    break;
		}
		case 5:
		{
		    color = "magenta";
		    break;
		}
		case 6:
		{
		    color = "yellow";
		    break;
		}
		case 7:
		{
		    color = "black";
		    break;
		}
		case 8:
		{
		    color = "indigo";
		    break;
		}
		case 9:
		{
		    color = "olive";
		    break;
		}
		case 10:
		{
		    color = "aqua";
		    break;
		}
		case 11:
		{
		    color = "violet";
		    break;
		}			    	    	    
		default:
		{
		    break;
		}
		}
		
		preview_svg.draw(display_polygon, color);
	    }
	    
	    preview_svg.Close();
	}
	else
	{
	    printf("Polygon optimization FAILED.\n");
	}
	getchar();
	
	vector<Polygon> next_polygons;
	
	for (unsigned int i = 0; i < remaining_polygons.size(); ++i)
	{
	    next_polygons.push_back(polygons[remaining_polygons[i]]);
	}
	
	polygon_index_map = remaining_polygons;
	polygons.clear();
	polygons = next_polygons;
    }
    while (!remaining_polygons.empty());

    finish = clock();
    
    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 15 ... finished\n");    
}    


void test_polygon_16(void)
{ 
    clock_t start, finish;
    
    printf("Testing polygon 16 ...\n");

    start = clock();

    SolverConfiguration solver_configuration;

    vector<Polygon> polygons;
    
    polygons.push_back(polygon_1);
    polygons.push_back(polygon_2);
    polygons.push_back(polygon_3);
    polygons.push_back(polygon_4);

    double area = calc_PolygonUnreachableZoneArea(polygon_1, polygons);
    printf("Polygons area: %.3f\n", area);

    finish = clock();
    
    printf("Time: %.3f\n", (finish - start) / (double)CLOCKS_PER_SEC);
    printf("Testing polygon 16 ... finished\n");    
}    


/*----------------------------------------------------------------*/

int main(int SEQ_UNUSED(argc), char **SEQ_UNUSED(argv))
{
    //test_polygon_1();
    //test_polygon_2();
    //test_polygon_3();
    //test_polygon_4();
    //test_polygon_5();
    //test_polygon_6();
    //test_polygon_7();
    //test_polygon_8();
    //test_polygon_9();
    //test_polygon_10();
    //test_polygon_11();
    //test_polygon_12();
    //test_polygon_13();
    //test_polygon_14();
    //test_polygon_15();
    test_polygon_16();    
    
    return 0;
}

