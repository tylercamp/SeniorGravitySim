


/*** Libraries ***/


#include <amp.h>
#include <amp_graphics.h>
#include <amp_math.h>
#include <vector>
#include <random>
#include <iostream>
#include <cassert>
#include <ShellScalingApi.h>


#include <Windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <SDL.h>

#include "Utilities.h"
#include "easybmp/EasyBMP.h"

#define NOT_YET_IMPLEMENTED() __debugbreak()

#pragma comment( lib, "opengl32.lib" )
#pragma comment( lib, "glu32.lib" )
#pragma comment( lib, "SDL2.lib" )



using concurrency::extent;

using namespace std;
using namespace concurrency;
using namespace fast_math;



#define RUN_BENCHMARK false
#define RUN_PRECISION_TEST false

#define DAYS(x) ((x) * 86400)
#define MONTHS(x) (DAYS(x * 30.42))
#define YEARS(x) (MONTHS(x * 12))

//	Bad name
#define FLOPS_PER_N(n) ((n) * (n) * 20 + (n) * 9)
#define GFLOPS(flops) ((flops) / 1e9)


#define TILE_SIZE 64

//#define WINDOW_WIDTH 6000
//#define WINDOW_HEIGHT 2000
//#define WINDOW_WIDTH (1920*1.8)
//#define WINDOW_HEIGHT (1080*1.8)
#define WINDOW_WIDTH ((int)(1680*1.5))
#define WINDOW_HEIGHT ((int)(1050*1.5))

#define USE_MASSIVE_SIM true
#define USING_NVIDIA true
#define USING_APPROXIMATION false
#define USING_FINE_STEP true

#define USE_NVIDIA_OPT true

#ifdef NDEBUG
//	Determine sim size from hardware params (don't want to freeze the GPU drivers)
#if USING_NVIDIA
# if USE_NVIDIA_OPT
#  if USE_MASSIVE_SIM
#   if USING_APPROXIMATION
#    define BODY_COUNT 30000
#   else
#    define BODY_COUNT 10000
#   endif
#  else
#   if USE_APPROXIMATION
#    define BODY_COUNT 15000
#   else
#    define BODY_COUNT 5000
#   endif
#  endif
# else
#  if USE_MASSIVE_SIM
#   define BODY_COUNT 5000
#  else
#   define BODY_COUNT 3000
#  endif
# endif
#else
# if USE_MASSIVE_SIM
#  if USING_APPROXIMATION
#   define BODY_COUNT 50000
#  else
#   define BODY_COUNT 5000
#  endif
# else
#  define BODY_COUNT 500
# endif
#endif
#else
#define BODY_COUNT 500
#endif

#if RUN_PRECISION_TEST
# define BODY_COUNT 5000
#endif


/* SIM PARAMS */
#ifndef EXPORT_DLL
# define SPATIAL_SCALE 1
#else
# define SPATIAL_SCALE 1e2
#endif

#define INITIAL_AREA 300

#define MOVE_SPEED 0.1f

//#define INITIAL_ENERGY 1e12
//#define INITIAL_ENERGY (BODY_COUNT*1e-3)
#define INITIAL_ENERGY 0

enum StartShape
{
	HollowSphere,
	MostlyHollowSphere,
	FilledSphere
};

#define INIT_HOLLOW_SPHERE 1
#define INIT_MOSTLY_HOLLOW_SPHERE 2
#define INIT_FILLED_SPHERE 3
#define INIT_FILLED_CUBE 4

#define INIT_TYPE INIT_HOLLOW_SPHERE





#define DEGTORAD(x) (x * 3.14159265358979323 / 180.0)

#ifdef EXPORT_DLL
# define API_FUNCTION extern "C" __declspec(dllexport)
#else
# define API_FUNCTION
#endif



   
#define get(type, name) __declspec(property(get=get_##name)) type name; \
	inline type get_##name( )

/*
struct body_data
{
	FLOAT radius;
	FLOAT mass;

	FLOAT_3 pos;
	FLOAT_3 vel;
	FLOAT_3 acc;

	body_data( )
	{
		radius = pow( ufrand( ), FLOAT( 4.0 ) ) * FLOAT( 1000.0 ) + FLOAT( 10.0 );
		mass = 4.0 * pow( radius, 3.0 ); //	approx sphere volume

		pos = FLOAT_3( frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA );
		vel = frand3( ) * sqrt( 2.0 * INITIAL_ENERGY / mass / BODY_COUNT );
	}
};
*/


enum TimeScale
{
	Seconds = 0,
	Days = 1,
	Months = 2,
	Years = 3
};

FLOAT convertTime( FLOAT time, TimeScale timeScale )
{
	switch( timeScale )
	{
	case(Seconds) :
		return time;
	case(Days) :
		return DAYS( time );
	case(Months) :
		return MONTHS( time );
	case(Years) :
		return YEARS( time );
	}
}

struct SimulationContext
{
	static SimulationContext * NewCopy( const SimulationContext & ref )
	{
		auto newCopy = new SimulationContext( );
		newCopy->timeScale = ref.timeScale;
		newCopy->timeStep = ref.timeStep;
		newCopy->body_positions = bindless_copy( ref.body_positions );
		newCopy->body_velocities = bindless_copy( ref.body_velocities );
		newCopy->body_accels = bindless_copy( ref.body_accels );

		return newCopy;
	}

	SimulationContext( ) :
		timeStep( 5e4 ),
		timeScale( Years ),
		body_positions( nullptr ),
		body_velocities( nullptr ),
		body_accels( nullptr ),
		GRAVITATIONAL_CONSTANT( FLOAT( 6.67e-11 ) )
	{
	}

	~SimulationContext( )
	{
		if( body_positions ) delete body_positions;
		if( body_velocities ) delete body_velocities;
		if( body_accels ) delete body_accels;
	}

	FLOAT GRAVITATIONAL_CONSTANT;

	FLOAT timeStep;
	TimeScale timeScale;

	array_view<FLOAT_4> * body_positions;
	array_view<FLOAT_3> * body_velocities;
	array_view<FLOAT_3> * body_accels;

	get( int, numBodies )
	{
		return body_positions->extent.size( );
	}
};



FLOAT_3 min3( const FLOAT_3 & a, const FLOAT_3 & b )
{
	return FLOAT_3(
		min( a.x, b.x ),
		min( a.y, b.y ),
		min( a.z, b.z )
		);
}

FLOAT_3 max3( const FLOAT_3 & a, const FLOAT_3 & b )
{
	return FLOAT_3(
		max( a.x, b.x ),
		max( a.y, b.y ),
		max( a.z, b.z )
		);
}




API_FUNCTION SimulationContext * start_new( )
{
	SimulationContext * result = new SimulationContext( );
	return result;
}

FLOAT_3 center_of_system( array_view<FLOAT_4> & bodies )
{
	int numBodies = bodies.extent.size( );

	array_view<FLOAT> body_localities( numBodies );
	body_localities.discard_data( );

	concurrency::extent<2> work_extent( numBodies, numBodies );

	parallel_for_each(
		work_extent,
		[=]( index<2> idx ) restrict( cpu, amp )
	{
		if( idx[0] == idx[1] )
			return;

		FLOAT_4 us = bodies[idx[0]];
		FLOAT_4 them = bodies[idx[1]];

		body_localities[idx[0]] += FLOAT( 1.0 ) / mag2( us.xyz - them.xyz );
	} );

	body_localities.synchronize( );


	
	FLOAT_3 center = FLOAT_3( 0.0 );
	FLOAT sumLocality = FLOAT( 0.0 );
	/*
	parallel_for_each(
		accelerator( accelerator::cpu_accelerator ).default_view,
		body_localities.extent.tile<128>( ).pad( ),
		[=, &center, &sumLocality]( tiled_index<1> idx ) restrict( cpu, amp )
	{
		if( idx.global[0] >= numBodies )
			return;

		FLOAT_3 tile_sum_position( 0.0 );
		FLOAT tile_sum_locality( 0.0 );

		for( int i = 0; i < 128; i++ )
		{
			int bodyIndex = idx.tile[0] * 128 + i;
			tile_sum_position += bodies[bodyIndex].pos;
			tile_sum_locality += body_localities[bodyIndex];
		}

		center += tile_sum_position;
		sumLocality += tile_sum_locality;
	} );
	*/

	for( int i = 0; i < numBodies; i++ )
	{
		center += bodies[i].xyz * body_localities[i];
		sumLocality += body_localities[i];
	}

	return center / sumLocality;
}

FLOAT_3 center_of_mass( SimulationContext * context )
{
	auto & bodies = *context->body_positions;

	FLOAT_3 netDisp = FLOAT_3( 0 );
	FLOAT sumMass = FLOAT( 0 );

	vector<FLOAT_4> data( bodies.extent.size( ) );
	for( int i = 0; i < bodies.extent.size( ); i++ )
		data[i] = bodies[i];

	int reduction_size = 10000;
	while( data.size( ) > reduction_size )
	{
		NOT_YET_IMPLEMENTED( );
	}

	for( int i = 0; i < data.size( ); i++ )
	{
		netDisp += data[i].xyz * data[i].w;
		sumMass += data[i].w;
	}

	return netDisp / sumMass;
}

API_FUNCTION FLOAT_3 range( SimulationContext * context )
{
	auto & bodies = *context->body_positions;

	FLOAT_3 min = FLOAT_3( 1e10 );
	FLOAT_3 max = FLOAT_3( -1e10 );

	vector<FLOAT_3> positions( bodies.extent.size( ) );
	for( int i = 0; i < bodies.extent.size( ); i++ )
		positions[i] = bodies[i].xyz;

	int reduction_size = 256;
	while( positions.size( ) > reduction_size )
	{
		NOT_YET_IMPLEMENTED( );
	}

	for( auto & pos : positions )
	{
		min = min3( pos, min );
		max = max3( pos, max );
	}

	return max - min;
}



API_FUNCTION void generate_bodies( SimulationContext * context, int count )
{
	if( context->body_positions )
		delete context->body_positions;
	if( context->body_velocities )
		delete context->body_velocities;
	if( context->body_accels )
		delete context->body_accels;

	auto & positions = context->body_positions;
	auto & velocities = context->body_velocities;
	auto & accels = context->body_accels;

	positions = new array_view<FLOAT_4>( count );
	velocities = new array_view<FLOAT_3>( count );
	accels = new array_view<FLOAT_3>( count );
	positions->discard_data( );
	velocities->discard_data( );
	accels->discard_data( );

	int numInnerBodies = count * 3 / 10;
	
	for( int i = 0; i < count; i++ )
	{
		FLOAT radius = 200.0;
		//FLOAT radius = pow( ufrand( ), FLOAT( 5.0 ) ) * FLOAT( 500.0 ) + FLOAT( 50.0 );
		FLOAT volume = 4.0 * pow( radius * SPATIAL_SCALE, 3.0 ); //	approx sphere volume
		FLOAT mass = volume;
		//FLOAT_3 pos = FLOAT_3( frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA );
		FLOAT_3 vel = frand3( ) * sqrt( 2.0 * INITIAL_ENERGY / mass / BODY_COUNT );

		FLOAT_3 pos;
#if INIT_TYPE == INIT_MOSTLY_HOLLOW_SPHERE
		if( i < numInnerBodies )
			pos = 0.1f * FLOAT_3( frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA );
		else
			pos = from_spherical( FLOAT_3( frand( ) * M_PI / 2, frand( ) * M_PI, INITIAL_AREA ) );
#endif
#if INIT_TYPE == INIT_HOLLOW_SPHERE
		pos = from_spherical( FLOAT_3( frand( ) * M_PI / 2, frand( ) * M_PI, INITIAL_AREA ) );
#endif
#if INIT_TYPE == INIT_FILLED_SPHERE
		pos = from_spherical( FLOAT_3( frand( ) * M_PI / 2, frand( ) * M_PI, INITIAL_AREA * sqrt( ufrand( ) ) ) );
#endif
#if INIT_TYPE == INIT_FILLED_CUBE
		pos = FLOAT_3( frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA, frand( ) * INITIAL_AREA );
#endif
#if INIT_TYPE == INIT_POINT_MATRIX
		int ix = (int)(ufrand( ) * ptMatrixDims);
		int iy = (int)(ufrand( ) * ptMatrixDims);
		int iz = (int)(ufrand( ) * ptMatrixDims);
		pos = FLOAT_3(
			ptMatrixStart + ix * ptMatrixStep + frand( ) * ptMatrixPtDims,
			ptMatrixStart + iy * ptMatrixStep + frand( ) * ptMatrixPtDims,
			ptMatrixStart + iz * ptMatrixStep + frand( ) * ptMatrixPtDims
		);
#endif

		positions->data( )[i] = FLOAT_4( pos.x, pos.y, pos.z, mass );
		velocities->data( )[i] = vel;
		accels->data( )[i] = FLOAT_3( 0.0 );
	}
}


API_FUNCTION void generate_bodies_uniform( SimulationContext * context, int count )
{
	if( context->body_positions )
		delete context->body_positions;
	if( context->body_velocities )
		delete context->body_velocities;
	if( context->body_accels )
		delete context->body_accels;

	auto & positions = context->body_positions;
	auto & velocities = context->body_velocities;
	auto & accels = context->body_accels;

	positions = new array_view<FLOAT_4>( count );
	velocities = new array_view<FLOAT_3>( count );
	accels = new array_view<FLOAT_3>( count );
	positions->discard_data( );
	velocities->discard_data( );
	accels->discard_data( );

	for( int i = 0; i < count; i++ )
	{
		FLOAT radius = 200.0;
		FLOAT volume = 4.0 * pow( radius * SPATIAL_SCALE, 3.0 ); //	approx sphere volume
		FLOAT mass = volume;
		FLOAT_3 vel = frand3_stable( ) * sqrt( 2.0 * INITIAL_ENERGY / mass / BODY_COUNT );

		FLOAT_3 pos = frand3_stable( ) * INITIAL_AREA;

		positions->data( )[i] = FLOAT_4( pos.x, pos.y, pos.z, mass );
		velocities->data( )[i] = vel;
		accels->data( )[i] = FLOAT_3( 0.0 );
	}
}

void draw_bodies( array_view<FLOAT_4> * bodies, FLOAT_3 centerOfMass, FLOAT_3 cameraPos )
{
	SDL_PumpEvents( );
	SDL_Event e;
	while( SDL_PollEvent( &e ) )
	{
		if( e.type == SDL_KEYDOWN )
		{
			//if( e.key )
		}
	}

	auto data = bodies->data( );
	
	glViewport( 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT );


	glEnable( GL_BLEND );
	glBlendFunc( GL_ONE, GL_ONE_MINUS_SRC_ALPHA );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );

	glBegin( GL_QUADS );
	//glColor4f( 0.0f, 0.0f, 0.0f, 1e-2f );
	glColor4f( 0.0f, 0.0f, 0.0f, 0.5f );
	glVertex2f( -1.0f, 1.0f );
	glVertex2f( 1.0f, 1.0f );
	glVertex2f( 1.0f, -1.0f );
	glVertex2f( -1.0f, -1.0f );
	glEnd( );

	glClear(GL_COLOR_BUFFER_BIT);

	glBlendFunc( GL_ONE_MINUS_DST_COLOR, GL_ONE );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	gluPerspective( 90.0f, WINDOW_WIDTH / (float)WINDOW_HEIGHT, 1e1, 1e20 );
	
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
	gluLookAt(
		-cameraPos.x, -cameraPos.y, -cameraPos.z,
		centerOfMass.x, centerOfMass.y, centerOfMass.z,
		0.0, 1.0, 0.0
	);
	glPointSize(1.5f);
	for( int i = 0; i < bodies->extent.size( ); i++ )
	{
		float b = 0.5f;
		glBegin( GL_POINTS );
		glColor3f( b, b, b );
		glVertex3d( data[i].x, data[i].y, data[i].z );
		glEnd( );
	}
}

API_FUNCTION void get_body_positions( SimulationContext * ctx, float * target, int count )
{
	FLOAT_4 * positions = ctx->body_positions->data( );
	for( int i = 0; i < count; i++ )
	{
		target[i * 3 + 0] = positions[i].x;
		target[i * 3 + 1] = positions[i].y;
		target[i * 3 + 2] = positions[i].z;
	}
}

API_FUNCTION void get_body_radii( SimulationContext * ctx, float * target, int count )
{
	FLOAT_4 * positions = ctx->body_positions->data( );
	for( int i = 0; i < count; i++ )
	{
		target[i] = pow( positions[i].w, 0.33333333 ) / 4.0;
	}
}

API_FUNCTION void set_body_radii( SimulationContext * ctx, float * source, int count )
{
	FLOAT_4 * bodies = ctx->body_positions->data( );
	for( int i = 0; i < count; i++ )
		bodies[i].w = 4.0 * pow( source[i], 3.0 );
}

API_FUNCTION void set_sim_timestep( SimulationContext * ctx, float timestep )
{
	ctx->timeStep = timestep;
}

API_FUNCTION void set_sim_time_scale( SimulationContext * ctx, int iTimeScale )
{
	TimeScale timeScale = (TimeScale)iTimeScale;
	ctx->timeScale = timeScale;
}

API_FUNCTION void set_grav_constant( SimulationContext * ctx, float constant )
{
	ctx->GRAVITATIONAL_CONSTANT = constant;
}




//	sampling optimization only works in this case due to complete randomization between spatial position and buffer index.
API_FUNCTION void simulate( SimulationContext * context, int samplePortions, int sampleIndex )
{
	auto & bodies = *context->body_positions;
	auto & vels = *context->body_velocities;
	auto & accels = *context->body_accels;

	auto numBodies = bodies.extent.size( );
	auto GRAVITATIONAL_CONSTANT = context->GRAVITATIONAL_CONSTANT;

#define EPS2 0.5

#if !USE_NVIDIA_OPT

	concurrency::extent<2> work_extent( numBodies, numBodies );

	parallel_for_each(
		work_extent,
		[=]( index<2> idx ) restrict( amp )
	{
		using namespace concurrency::fast_math;

		if( idx[0] == idx[1] )
			return;

		const FLOAT_4 a = bodies[idx[0]];
		const FLOAT_4 b = bodies[idx[1]];
		
		FLOAT_3 ab = b.xyz - a.xyz;
		ab *= SPATIAL_SCALE;
		FLOAT dist2 = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z + EPS2;
		FLOAT dist6 = dist2 * dist2 * dist2;
		FLOAT invDist3 = FLOAT( 1.0 ) / sqrt( dist6 );
		//float force = GRAVITATIONAL_CONSTANT * a.w * b.w / ( dist * dist + 1.0 );

		accels[idx[0]] += b.w * ab * invDist3;
		
	} );
	
#else

	if( numBodies % samplePortions != 0 )
		NOT_YET_IMPLEMENTED( );

	concurrency::extent<2> work_extent( numBodies, ceil( numBodies / (float) TILE_SIZE / samplePortions ) );

	int offset = numBodies / samplePortions * sampleIndex;

	parallel_for_each(
		work_extent.tile<TILE_SIZE, 1>().pad(),
		[=]( tiled_index<TILE_SIZE, 1> t_idx ) restrict( amp )
	{
		tile_static FLOAT_4 tile_bodies[TILE_SIZE];
		tile_bodies[t_idx.local[0]] = bodies[offset + t_idx.global[1] * TILE_SIZE + t_idx.local[0]];

		FLOAT_4 thread_body = bodies[t_idx.global[0]];
		t_idx.barrier.wait( );

		FLOAT_3 acc = FLOAT_3( 0.0 );

		for( int i = 0; i < TILE_SIZE; i++ )
		{
			if( i + t_idx.global[1] * TILE_SIZE >= numBodies )
				continue;

			const FLOAT_4 a = thread_body;
			const FLOAT_4 b = tile_bodies[i];

			FLOAT_3 ab = b.xyz - a.xyz;
			ab *= SPATIAL_SCALE;
			FLOAT dist2 = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z + EPS2;
			FLOAT dist6 = dist2 * dist2 * dist2;
			FLOAT invDist6 = FLOAT( 1.0 ) / dist6;
			FLOAT invDist3 = sqrt( invDist6 );

			acc += b.w * ab * invDist3 * GRAVITATIONAL_CONSTANT * samplePortions;
		}

		t_idx.barrier.wait( );
		accels[t_idx.global[0]] += acc;
	} );
#endif

	FLOAT step = convertTime( context->timeStep, context->timeScale );
	parallel_for_each(
		bodies.extent,
		[=]( index<1> idx ) restrict( amp )
	{
		//	TODO: Use some 4th-order approx.
		//	http://gafferongames.com/game-physics/integration-basics/

		FLOAT_4 b = bodies[idx];
		FLOAT_3 b_vel = vels[idx];
		FLOAT_3 b_acc = GRAVITATIONAL_CONSTANT * accels[idx] / ( b.w * numBodies ) ;
		b.xyz += 0.98 * b_vel * step;
		b_vel += 0.98 * b_acc * step;

		bodies[idx] = b;
		vels[idx] = b_vel;
		accels[idx] = FLOAT_3( 0.0 );
	} );
}

//	Precision compared to complete approximation for current method
void run_precision_test( SimulationContext * ctx )
{
	const auto len2 = []( const FLOAT_3 & v )
	{
		return v.x * v.x + v.y * v.y + v.z * v.z;
	};
	const auto sqrt = []( const FLOAT_3 & v )
	{
		return FLOAT_3( ::sqrt( v.x ), ::sqrt( v.y ), ::sqrt( v.z ) );
	};

	auto preciseSim = SimulationContext::NewCopy( *ctx );
	auto lossySim = SimulationContext::NewCopy( *ctx );
	auto & precisePositions = *preciseSim->body_positions;
	auto & lossyPositions = *lossySim->body_positions;

	int numSteps = 7500;
	int approximationConstant = 100;

	auto calculateSumError = [=]( )
	{
		double error = 0.0;
		for( int i = 0; i < ctx->numBodies; i++ )
		{
			FLOAT_4 diff = precisePositions[i] - lossyPositions[i];
			error += mag( diff.xyz );
		}
		return error;
	};
	auto calculatePreciseCenter = [=]( )
	{
		FLOAT_4 sumRealDist( 0.0 );
		for( int i = 0; i < ctx->numBodies; i++ )
			sumRealDist += precisePositions[i];
		return sumRealDist.xyz / ctx->numBodies;
	};
	auto calculatePreciseStdDev = [=]( FLOAT_3 mean )
	{
		FLOAT_3 variance = 0;
		for( int i = 0; i < ctx->numBodies; i++ )
			variance += (precisePositions[i].xyz - mean) * (precisePositions[i].xyz - mean);
		return sqrt( variance );
	};

	std::cout << "Generating up to " << numSteps << " steps" << std::endl;
	std::cout << "Approximating on the order of " << approximationConstant << std::endl;
	/*
	std::cout << "Generating precise reference... ";
	for( int i = 0; i < numSteps; i++ )
		simulate( preciseSim, 1, 0 );

	preciseSim->body_positions->synchronize( );
	std::cout << "Done." << std::endl;
	*/

	std::vector<float> imprecisionRecord;

	std::cout << "Running approximation...\n";
	for( int i = 0; i < numSteps; i++ )
	{
		simulate( lossySim, approximationConstant, i % approximationConstant );
		simulate( preciseSim, 1, i % 1 );

		//	Occasionally sync for responsiveness
		if( i % 5 == 0 )
		{
			auto acc = accelerator( accelerator::default_accelerator ).default_view;
			acc.flush( );
			acc.wait( );
			Sleep( 10 );
		}

		if( i % ( numSteps / 100 ) == 0 )
		{
			std::cout << i << ": ";
			FLOAT_3 stdDev = calculatePreciseStdDev( calculatePreciseCenter( ) );
			auto me_sd = calculateSumError( ) / ctx->numBodies / sqrtf( len2( stdDev ) );
			//std::cout << "(Mean error) / len(Std Dev.) = ";
			std::cout << "ME/SD = ";
			std::cout << me_sd << std::endl;
			imprecisionRecord.push_back( me_sd );
		}
	}
	std::cout << "Done." << std::endl;

	auto secondOrderErrorTrend = [=]( )
	{
		//	Sample 1
		int s1_a = 0, s1_b = 1;
		//	Sample 2
		int s2_a = imprecisionRecord.size( ) - 2, s2_b = imprecisionRecord.size( ) - 1;
		float ds1 = ( imprecisionRecord[s1_b] - imprecisionRecord[s1_a] ) / ( s1_b - s1_a );
		float ds2 = ( imprecisionRecord[s2_b] - imprecisionRecord[s2_a] ) / ( s2_b - s1_a );
		return ( ds2 - ds1 ) / ( s2_b - s1_b );
	};

	std::cout << "Simulated time: " << convertTime( numSteps * ctx->timeStep, ctx->timeScale )  / YEARS( 1 ) << " years" << std::endl;

	auto stdDev = calculatePreciseStdDev( calculatePreciseCenter( ) );
	printf( "Mean spatial dist from 0: %f %f %f\n", stdDev.x, stdDev.y, stdDev.z );
	std::cout << "Mean spatial error: " << calculateSumError( ) / ctx->numBodies  << std::endl;
	std::cout << "Error trend approximation: diff(error)/steps = " << (imprecisionRecord.back( ) - imprecisionRecord.front( )) / numSteps << std::endl;
	std::cout << "dd(error)/dds = " << secondOrderErrorTrend( ) << std::endl;
}

void run_benchmark( SimulationContext * ctx )
{
	auto accel_view = accelerator( accelerator::default_accelerator ).default_view;
	auto start = clock( );
	int numSteps = 100;
	for( int i = 0; i < numSteps; i++ )
		simulate( ctx, 1, 0 );

	accel_view.wait( );
	auto end = clock( );
	auto seconds = (end - start) / 1000.0;
	std::cout << "GFLOPS: " << GFLOPS( FLOPS_PER_N( BODY_COUNT ) ) * numSteps / seconds << std::endl;
}

void save_screenshot(const std::string & file, std::uint8_t * screen)
{
	BMP bmp;
	bmp.SetSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	bmp.SetBitDepth(24);

	for(int y = 0; y < WINDOW_HEIGHT; y++)
	{
		for(int x = 0; x < WINDOW_WIDTH; x++)
		{
			auto start = screen + 3 * (y * WINDOW_WIDTH + x);
			bmp.SetPixel(x, y, { start[0], start[1], start[2], 255 });
		}
	}
	bmp.WriteToFile(file.c_str());
}

#ifndef DLL_EXPORT
#undef main
int main( )
{
	SetProcessDPIAware( );

	g_ThreadPool = new ThreadPool( 4 );

	auto context = start_new( );

	generate_bodies( context, BODY_COUNT );

#if RUN_BENCHMARK
	run_benchmark( context );
	system( "pause" );

	delete g_ThreadPool;
	delete context;
	return 0;
#endif

#if RUN_PRECISION_TEST
	generate_bodies_uniform( context, BODY_COUNT );
	run_precision_test( context );
	system( "pause" );

	delete g_ThreadPool;
	delete context;
	return 0;
#endif

	SDL_Init( SDL_INIT_VIDEO );

	SDL_Window * window = SDL_CreateWindow( "Gravity", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_OPENGL );
	SDL_GL_CreateContext( window );

	//	Spherical coordinates
	static FLOAT_3 offset = FLOAT_3( DEGTORAD( -45.0 ), DEGTORAD( -45.0 ), 1e4 );

	auto t = clock( );
	int frames = 0;
	int tick = 0;

	//	1 = exact approximation
	const int approximationConstant = 250;

	std::size_t screenshotIdx = 0;
	std::uint8_t * screenBuffer = new std::uint8_t[3 * WINDOW_WIDTH * WINDOW_HEIGHT];


	FLOAT numSeconds = 0.0;
	int prevMagnitude = 0;
	bool run = true;
	while( run )
	{
		if( clock( ) - t > 1000 )
		{
			float approximationTime = approximationConstant / (float)frames;
			printf( "%ifps (%f yrs per second, approximation time = %fs)\n", frames, numSeconds / YEARS(1), approximationTime );
			t = clock( );
			frames = 0;
			numSeconds = FLOAT( 0.0 );
		}

		++frames;
		++tick;

		FLOAT_3 prev_offset = offset;
		if( !GetAsyncKeyState( VK_SPACE ) )
			offset.x += 0.001;

		if( GetAsyncKeyState( VK_LEFT ) ) offset.x -= 0.01;
		if( GetAsyncKeyState( VK_RIGHT ) ) offset.x += 0.01;
		if( GetAsyncKeyState( VK_UP ) ) offset.z += 0.02 * offset.z;
		if( GetAsyncKeyState( VK_DOWN ) ) offset.z -= 0.02 * offset.z;

		if( prev_offset != offset )
			glClear( GL_COLOR_BUFFER_BIT );

		//	Left brace [
		if( GetAsyncKeyState( 219 ) ) context->timeStep *= FLOAT( 0.95 );
		//	Right brace ]
		if( GetAsyncKeyState( 221 ) ) context->timeStep /= FLOAT( 0.95 );

		if( GetAsyncKeyState( VK_ESCAPE ) )
			run = false;

		//FLOAT_3 centerOfMass = center_of_system( points );
		//FLOAT_3 centerOfMass = center_of_mass( context );
		//if( GetAsyncKeyState( VK_SPACE ) ) 
		//	centerOfMass = FLOAT_3( 0.0 );
		//FLOAT_3 cameraPosition = centerOfMass + from_spherical( offset );
		FLOAT_3 cameraPosition = from_spherical( offset );
		int magnitude = static_cast<int>(log10( offset.z * SPATIAL_SCALE ));
		if( prevMagnitude != magnitude )
		{
			std::cout << "Magnitude: 1 * 10^" << magnitude << "\n";
		}

		prevMagnitude = magnitude;


#if USING_FINE_STEP
		context->timeStep /= 1000;
		for(int i = 0; i < 1000; i++)
		{
#endif
#if USING_APPROXIMATION
			simulate(context, approximationConstant, tick % approximationConstant);
#else
			simulate(context, 1, tick % 1);
#endif
#if USING_FINE_STEP
			Sleep(1);
		}
		context->timeStep *= 1000;
#endif
		numSeconds += convertTime( context->timeStep, context->timeScale );


		draw_bodies( context->body_positions, FLOAT_3( 0.0 ), cameraPosition );
		SDL_GL_SwapWindow( window );

#if USING_FINE_STEP
		glReadPixels(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, screenBuffer);
		save_screenshot("frames/frame_" + tostring(++screenshotIdx) + ".bmp", screenBuffer);
#endif


		//Sleep( 50 );
		Sleep( 1 );
	}

	SDL_Quit( );

	delete g_ThreadPool;
	delete context;

	//std::cout << "Total frames: " << tick << std::endl;
	//system( "pause" );

	return 0;
}
#endif