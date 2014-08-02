/* File: timer.h
 *
 * Provides an abstraction for a Timer object to keep track of the processing
 * time of the various algorithm implementation. This class allows different
 * libraries to be swapped out in the underlying timer implementation.
 */

#ifndef TIMER_H
#define TIMER_H

#include "boost/date_time/local_time/local_time.hpp"


class Timer {
  private:
    boost::posix_time::ptime start_t;
    bool started; // true if the timer was ever started with start()
    bool running; // true if the timer is currently ticking
    unsigned long counter;

    // returns the number of milliseconds since start() was called
    unsigned long getDiff();


  public:
    Timer();

    // start and stop the clock
    void start();
    void stop();

    // reset the counter to 0
    void reset();

    // return the number of milliseconds that passed
    unsigned long get();
};


#endif
