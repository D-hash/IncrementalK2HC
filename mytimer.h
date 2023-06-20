#ifndef TIMER_INC
#define TIMER_INC
#include <chrono>
#include <sys/time.h>
#include <sys/resource.h> 
#include <unistd.h> 
#include <boost/timer/timer.hpp>
#include <cmath>

class mytimer{
//Redefining to avoid conflicts with PartitionMETIS libraries (timer data type is a double)
public:
	mytimer() {
		start = timestamp(); // postcondition: elapsed()==0
		pauseSum=0;
  	} 
	void restart() { 
		start = timestamp(); // post: elapsed()==0
		pauseSum=0;
  	}
	void pause(){
		wpause=timestamp();
	}
	void resume(){
		pauseSum += (timestamp()- wpause);
	}
  	double elapsed(){                  // return elapsed time in seconds
		return timestamp()-start - pauseSum;
	}

private:
	double start;
	double wpause;
	double pauseSum;

	inline double timestamp() {
	
//	return std::chrono::high_resolution_clock::now();
	//Precise timing
	 struct timespec time;
	//CLOCK_PROCESS_CPUTIME_ID replaced with CLOCK_MONOTONIC to deal with multicore architecture

	clock_gettime(CLOCK_MONOTONIC, &time);
	 return (double)time.tv_sec + (double)time.tv_nsec * 1e-9;
//

	// Old timing
	//	struct timeval tp;
	//	gettimeofday(&tp, NULL);
	//	return double(tp.tv_sec) + tp.tv_usec / 1000000.;
	
	}
	
};




#endif  // TIMER_INC
