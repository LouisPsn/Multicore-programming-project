#ifndef PERFCOUNTER_IS_DEF
#define PERFCOUNTER_IS_DEF

#include <stdint.h>
#include <stdbool.h>

// Monitor flags
#define EASYPAP_MONITOR_L1 1
#define EASYPAP_MONITOR_L2 2
#define EASYPAP_MONITOR_L3 4
#define EASYPAP_MONITOR_TLB 8
#define EASYPAP_MONITOR_ALL                                                    \
  (EASYPAP_MONITOR_L1 | EASYPAP_MONITOR_L2 | EASYPAP_MONITOR_L3 |              \
   EASYPAP_MONITOR_TLB)

typedef enum
{
  EASYPAP_L2_HIT,
  EASYPAP_L3_MISS,
  EASYPAP_L3_HIT,
  EASYPAP_ALL_LOADS,
  EASYPAP_NB_COUNTERS
} easypap_perfcounter_counter_t;

typedef struct
{
  bool is_native;
  union {
    char *name;
    int code;
  };
} perf_event;

#define ERROR_RETURN(retval)                                                   \
  {                                                                            \
    fprintf (stderr, "Error %d %s:line %d: \n", retval, __FILE__, __LINE__);   \
    exit (retval);                                                             \
  }

void easypap_perfcounter_init (unsigned nb_cpus, unsigned monitor_flags);

void easypap_perfcounter_monitor_start_tile (unsigned cpu);
void easypap_perfcounter_monitor_stop_tile (unsigned cpu);

void easypap_perfcounter_monitor_start_iteration ();
void easypap_perfcounter_monitor_stop_iteration ();

int easypap_perfcounter_create_event_set (unsigned cpu);
int64_t easypap_perfcounter_get_counter (unsigned cpu,
                                         easypap_perfcounter_counter_t counter);
int easypap_perfcounter_get_counters (int64_t *counter_array, unsigned cpu);
int easypap_perfcounter_get_total_counters (int64_t *total_counters);

void easypap_perfcounter_finalize (void);

#if ENABLE_PAPI
extern unsigned do_cache;
#else
#define do_cache (unsigned)0
#endif

#endif
