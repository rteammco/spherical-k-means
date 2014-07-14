/* File: timer.h
 *
 * Provides an abstraction for a Timer object to keep track of the processing
 * time of the various algorithm implementation. This class allows different
 * libraries to be swapped out in the underlying timer implementation.
 */

#ifndef TIMER_H
#define TIMER_H

#include "Galois/Timer.h"


class Timer {
  private:
    Galois::Timer galois_timer;
    unsigned long counter;

  public:
    Timer();

    // start and stop the clock
    void start();
    void stop();

    // reset the counter to 0
    void reset();

    // return the number of miliseconds that passed
    unsigned long get();
};


#endif
