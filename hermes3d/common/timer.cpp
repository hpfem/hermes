// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

//
// timer.cc
//

/*
#include <config.h>
#include <common/memcheck.h>
#ifdef DEBUG
#define new DEBUG_NEW
#define malloc DEBUG_MALLOC
#define free DEBUG_FREE
#endif
*/
#include "timer.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// Timer //////////////////////////////////////////////////////////////////////

Timer::Timer() {
        tick(H3D_SKIP);
	reset();
	this->name = NULL;
}

Timer::Timer(const char *name) {
        tick(H3D_SKIP);
	reset();
	this->name = new char [strlen(name) + 1];
	strcpy(this->name, name);
}

Timer::~Timer() {
	delete [] name;
}

void Timer::start(bool rst) {
	if (rst) reset();
	gettimeofday(&this->st, NULL);
}

void Timer::tick(TimerTickType type) {
  if (type == H3D_SKIP) {
    start();
    accum_time = 0;
    last_period = -1;
  }
  else {
    stop(); 
    double secs = get_seconds();
    accum_time += secs;
    last_period = secs;
  }
}

void Timer::stop() {
	struct timeval tv;
	gettimeofday(&tv, NULL);

	accum.tv_sec += tv.tv_sec - st.tv_sec;
	accum.tv_usec += tv.tv_usec - st.tv_usec;
}

void Timer::reset() {
	timerclear(&this->accum);
}

double Timer::get_seconds() {
	return accum.tv_sec + (accum.tv_usec / 1000000.0);
}

const char *Timer::get_human_time() {
	#define MAX_STR					200
	// FIXME: not thread safe!
	static char str[MAX_STR] = "";
	double secs = get_seconds();
	int hours = (int) secs / (3600);
	int mins = (int) fmod(secs, 3600) / 60;
	sprintf(str, "%dh %dm %02.2lfs", hours, mins, fmod(fmod(secs, 3600), 60));

	return str;
}

