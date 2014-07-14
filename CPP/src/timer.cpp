/* File: timer.cpp
 *
 * Provides an abstraction for a Timer object to keep track of the processing
 * time of the various algorithm implementation. This class allows different
 * libraries to be swapped out in the underlying timer implementation.
 */

#include "timer.h"


// Constructor: initializes the timer's counter to 0.
Timer::Timer()
{
    counter = 0;
}


// Start the timer. The timer will start where it left off (or 0 initially
// or after a reset).
void Timer::start()
{
    galois_timer.start();
}


// Stops the timer from ticking. This does not reset the counter.
void Timer::stop()
{
    galois_timer.stop();
    counter += galois_timer.get();
}


// Resets the timer's counter back to 0.
void Timer::reset()
{
    counter = 0;
}


// Returns the number of seconds elapsed while the timer was running. Stopping
// the timer does not reset the counter value.
unsigned long Timer::get()
{
    return counter;
}
