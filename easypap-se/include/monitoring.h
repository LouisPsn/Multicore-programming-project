#ifndef MONITORING_IS_DEF
#define MONITORING_IS_DEF

#include "gmonitor.h"
#include "perfcounter.h"
#include "time_macros.h"
#include "trace_record.h"

#ifdef ENABLE_MONITORING

static inline void monitoring_declare_task_ids (char *task_ids[])
{
  trace_record_declare_task_ids (task_ids);
}

#ifdef ENABLE_SDL

static inline long monitoring_start_iteration (void)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_start_iteration ();
#endif
  if (do_gmonitor | do_trace) {
    long t = what_time_is_it ();
    gmonitor_start_iteration (t);
    trace_record_start_iteration (t);
    return t;
  } else
    return 0;
}

static inline long monitoring_end_iteration (void)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_stop_iteration ();
#endif
  if (do_gmonitor | do_trace) {
    long t = what_time_is_it ();
    gmonitor_end_iteration (t);
    trace_record_end_iteration (t);
    return t;
  } else
    return 0;
}

static inline void monitoring_start_tile (unsigned cpu)
{
  if (do_gmonitor | do_trace) {
    long t = what_time_is_it ();
    gmonitor_start_tile (t, cpu);
    trace_record_start_tile (t, cpu);
  }
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_start_tile (cpu);
#endif
}

static inline void monitoring_end_tile (unsigned x, unsigned y, unsigned w,
                                        unsigned h, unsigned cpu)
{
#ifdef ENABLE_PAPI

  if (do_cache)
    easypap_perfcounter_monitor_stop_tile (cpu);
#endif
  if (do_gmonitor | do_trace) {
    long t = what_time_is_it ();
    gmonitor_end_tile (t, cpu, x, y, w, h);
#ifdef ENABLE_PAPI
    if (do_cache) {
      int64_t counters[EASYPAP_NB_COUNTERS];
      easypap_perfcounter_get_counters (counters, cpu);
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0,
                             counters);
    } else {
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0, NULL);
    }
#else
    trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0, NULL);
#endif
  }
}

static inline void monitoring_end_tile_id (unsigned x, unsigned y, unsigned w,
                                           unsigned h, unsigned cpu,
                                           unsigned task_id)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_stop_tile (cpu);
#endif
  if (do_gmonitor | do_trace) {
    long t = what_time_is_it ();
    gmonitor_end_tile (t, cpu, x, y, w, h);
#ifdef ENABLE_PAPI
    if (do_cache) {
      int64_t counters[EASYPAP_NB_COUNTERS];
      easypap_perfcounter_get_counters (counters, cpu);
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                             counters);
    } else {

      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                             NULL);
    }
#else
    trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                           NULL);
#endif
  }
}

static inline void monitoring_gpu_tile (unsigned x, unsigned y, unsigned w,
                                        unsigned h, unsigned cpu, long start,
                                        long end, task_type_t task_type)
{
  if (do_gmonitor | do_trace) {
    if (task_type == TASK_TYPE_COMPUTE)
      gmonitor_start_tile (start, cpu);
    trace_record_start_tile (start, cpu);
    if (task_type == TASK_TYPE_COMPUTE)
      gmonitor_end_tile (end, cpu, x, y, w, h);
    trace_record_end_tile (end, cpu, x, y, w, h, task_type, 0, NULL);
  }
}

#else // no SDL

static inline long monitoring_start_iteration (void)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_start_iteration ();
#endif
  if (do_trace) {
    long t = what_time_is_it ();
    trace_record_start_iteration (t);
    return t;
  } else
    return 0;
}

static inline long monitoring_end_iteration (void)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_stop_iteration ();
#endif
  if (do_trace) {
    long t = what_time_is_it ();
    trace_record_end_iteration (t);
    return t;
  } else
    return 0;
}

static inline void monitoring_start_tile (unsigned cpu)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_start_tile (cpu);
#endif
  if (do_trace) {
    long t = what_time_is_it ();
    trace_record_start_tile (t, cpu);
  }
}

static inline void monitoring_end_tile (unsigned x, unsigned y, unsigned w,
                                        unsigned h, unsigned cpu)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_stop_tile (cpu);
#endif
  if (do_trace) {
    long t = what_time_is_it ();
#ifdef ENABLE_PAPI
    if (do_cache) {
      int64_t counters[EASYPAP_NB_COUNTERS];
      easypap_perfcounter_get_counters (counters, cpu);
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0,
                             counters);
    } else {
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0, NULL);
    }
#else
    trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, 0, NULL);
#endif
  }
}

static inline void monitoring_end_tile_id (unsigned x, unsigned y, unsigned w,
                                           unsigned h, unsigned cpu,
                                           unsigned task_id)
{
#ifdef ENABLE_PAPI
  if (do_cache)
    easypap_perfcounter_monitor_stop_tile (cpu);
#endif
  if (do_trace) {
    long t = what_time_is_it ();
#ifdef ENABLE_PAPI
    if (do_cache) {
      int64_t counters[EASYPAP_NB_COUNTERS];
      easypap_perfcounter_get_counters (counters, cpu);
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                             counters);
    } else {
      trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                             NULL);
    }
#else
    trace_record_end_tile (t, cpu, x, y, w, h, TASK_TYPE_COMPUTE, task_id + 1,
                           NULL);
#endif
  }
}

static inline void monitoring_gpu_tile (unsigned x, unsigned y, unsigned w,
                                        unsigned h, unsigned cpu, long start,
                                        long end, task_type_t task_type)
{
  if (do_trace) {
    trace_record_start_tile (start, cpu);
    trace_record_end_tile (end, cpu, x, y, w, h, task_type, 0, NULL);
  }
}

#endif

#else

#define monitoring_declare_task_ids (task_ids) (void) 0
#define monitoring_start_iteration() (void)0
#define monitoring_end_iteration() (void)0
#define monitoring_start_tile(c) (void)0
#define monitoring_end_tile(x, y, w, h, c) (void)0
#define monitoring_end_tile_id(x, y, w, h, c, id) (void)0
#define monitoring_gpu_tile(x, y, w, h, c, s, e, tt) (void)0

#endif

#endif
