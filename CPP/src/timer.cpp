/* File: timer.cpp
 *
 * Provides an abstraction for a Timer object to keep track of the processing
 * time of the various algorithm implementation. This class allows different
 * libraries to be swapped out in the underlying timer implementation.
 */

#include "timer.h"

using namespace boost::posix_time;


// Constructor: initializes the timer's counter to 0.
Timer::Timer()
{
    started = false;
    running = false;
    counter = 0;
}


// Start the timer. The timer will start where it left off (or 0 initially
// or after a reset).
void Timer::start()
{
    start_t = microsec_clock::local_time();
    started = true;
    running = true;
}


// Stops the timer from ticking. This does not reset the counter.
void Timer::stop()
{
    if(running) {
        counter += getDiff();
        running = false;
    }
}


// Resets the timer's counter back to 0.
void Timer::reset()
{
    started = false;
    counter = 0;
}


// Returns the number of milliseconds that passed since the last start()
// functon call. If it was never called, returns 0.
unsigned long Timer::getDiff()
{
    if(!started)
        return 0;

    ptime end_t = microsec_clock::local_time();
    time_duration diff = end_t - start_t;
    return diff.total_milliseconds();
}


// Returns the number of seconds elapsed while the timer was running. Stopping
// the timer does not reset the counter value.
unsigned long Timer::get()
{
    if(running)
        return getDiff();
    else
        return counter;
}
