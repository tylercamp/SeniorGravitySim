#pragma once

#include <amp.h>
#include <amp_graphics.h>
#include <amp_math.h>
#include <vector>
#include <random>
#include <iostream>


#include <Windows.h>

using namespace concurrency;

#define USE_HIGH_PRECISION true

//typedef float FLOAT;
//typedef graphics::double_3 FLOAT_3;

#if !USE_HIGH_PRECISION
# define FLOAT float
# define FLOAT_3 graphics::float_3
# define FLOAT_4 graphics::float_4
#else
# define FLOAT double
# define FLOAT_3 graphics::double_3
# define FLOAT_4 graphics::double_4
#endif



#pragma region Thread Pool

//	https://github.com/progschj/ThreadPool

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <random>

class ThreadPool
{
public:
	ThreadPool( size_t );
	template<class F, class... Args>
	auto Enqueue( F&& f, Args&&... args )
		->std::future<typename std::result_of<F( Args... )>::type>;
	void WaitForAll( );
	~ThreadPool( );
private:
	// need to keep track of threads so we can join them
	std::vector< std::thread > workers;
	// the task queue
	std::queue< std::function<void( )> > tasks;

	// synchronization
	std::mutex queue_mutex;
	std::condition_variable condition;
	bool stop;
	bool prepareForReset;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool( size_t threads )
	: stop( false ), prepareForReset( false )
{
	for( size_t i = 0; i<threads; ++i )
		workers.emplace_back(
			[this]
	{
		for( ;;)
		{
			std::function<void( )> task;

			{
				std::unique_lock<std::mutex> lock( this->queue_mutex );
				this->condition.wait( lock,
					[this] { return this->stop || !this->tasks.empty( ); } );
				if( this->tasks.empty( ) )
				{
					if( this->stop )
						return;
					if( this->prepareForReset )
					{
						while( this->prepareForReset )
							Sleep( 1 );
					}
				}
				task = std::move( this->tasks.front( ) );
				this->tasks.pop( );
			}

			task( );
		}
	}
	);
}

// add new work item to the pool
template<class F, class... Args>
auto ThreadPool::Enqueue( F&& f, Args&&... args )
-> std::future<typename std::result_of<F( Args... )>::type>
{
	using return_type = typename std::result_of<F( Args... )>::type;

	auto task = std::make_shared< std::packaged_task<return_type( )> >(
		std::bind( std::forward<F>( f ), std::forward<Args>( args )... )
		);

	std::future<return_type> res = task->get_future( );
	{
		std::unique_lock<std::mutex> lock( queue_mutex );

		// don't allow enqueueing after stopping the pool
		if( stop )
			throw std::runtime_error( "enqueue on stopped ThreadPool" );

		tasks.emplace( [task]( ) { (*task)(); } );
	}
	condition.notify_one( );
	return res;
}

inline void ThreadPool::WaitForAll( )
{

}

// the destructor joins all threads
inline ThreadPool::~ThreadPool( )
{
	{
		std::unique_lock<std::mutex> lock( queue_mutex );
		stop = true;
	}
	condition.notify_all( );
	for( std::thread &worker : workers )
		worker.join( );
}

#pragma endregion


ThreadPool * g_ThreadPool = nullptr;


/*** Utilities ***/

#undef max
#undef min

#define STABLE_RANDOM_SEED 1234

FLOAT ufrand( )
{
	static std::mt19937 random( clock( ) );
	return static_cast<double>(random( )) / random.max( ); // range 0 to 1
}

FLOAT ufrand_stable( )
{
	static std::mt19937 random_stable( STABLE_RANDOM_SEED );
	return static_cast<double>( random_stable( ) ) / random_stable.max( ); // range 0 to 1
}

FLOAT frand( )
{
	return ufrand( ) * FLOAT( 2.0 ) - FLOAT( 1.0 ); // range -1 to 1
}

FLOAT frand_stable( )
{
	return ufrand_stable( ) * FLOAT( 2.0 ) - FLOAT( 1.0 ); // range -1 to 1
}

//	x, y, z, as theta, inc, r
FLOAT_3 from_spherical( const FLOAT_3 & spherical )
{
	float rot = spherical.x;
	float inc = spherical.y;
	float r = spherical.z;
	return r * FLOAT_3(
		cos( rot ) * cos( inc ),
		sin( inc ),
		sin( rot ) * cos( inc )
		);
}

FLOAT_3 frand3( )
{
	FLOAT theta, inc;
	theta = frand( ) * FLOAT( 180.0 );
	inc = frand( ) * FLOAT( 90.0 );

	return from_spherical( FLOAT_3( theta, inc, 1.0 ) );
}

FLOAT_3 frand3_stable( )
{
	FLOAT theta, inc;
	theta = frand_stable( ) * FLOAT( 180.0 );
	inc = frand_stable( ) * FLOAT( 90.0 );

	return from_spherical( FLOAT_3( theta, inc, 1.0 ) );
}

int irand( )
{
	static std::mt19937 random( clock( ) );
	return random( );
}

int irand_stable( )
{
	static std::mt19937 random( STABLE_RANDOM_SEED );
	return random( );
}

FLOAT mag2( const FLOAT_3 & v ) restrict( cpu )
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

FLOAT mag2( const FLOAT_3 & v ) restrict( amp )
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

FLOAT mag( const FLOAT_3 & v ) restrict( cpu )
{
	return sqrt( mag2( v ) );
}

FLOAT mag( const FLOAT_3 & v ) restrict( amp )
{
	return fast_math::sqrt( mag2( v ) );
}

template <typename T>
std::string tostring(T num)
{
	std::ostringstream strm;
	strm << num;
	return strm.str();
}

template <typename T>
concurrency::array_view<T> * bindless_copy( const std::vector<T> & source )
{
	auto result = new concurrency::array_view<T>( (int) source.size( ) );
	result->discard_data( );

	for( std::size_t i = 0; i < source.size( ); i++ )
	{
		( *result )[i] = source[i];
	}

	return result;
}

template <typename T>
concurrency::array_view<T> * bindless_copy( const concurrency::array_view<T> & source )
{
	auto result = new concurrency::array_view<T>( source.extent );
	result->discard_data( );

	for( std::size_t i = 0; i < source.extent.size( ); i++ )
	{
		( *result )[i] = source[i];
	}

	return result;
}

template <typename T>
inline concurrency::array_view<T> * bindless_copy( const concurrency::array_view<T> const * source )
{
	return bindless_copy( *source );
}