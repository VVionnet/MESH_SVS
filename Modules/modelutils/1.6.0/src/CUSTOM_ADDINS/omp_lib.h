! default integer type assumed below
! default logical type assumed below
! OpenMP API v4.5

include ’omp_lib_kinds.h’
integer openmp_version
parameter ( openmp_version = 201511 )

external omp_set_num_threads
external omp_get_num_threads
integer omp_get_num_threads
external omp_get_max_threads
integer omp_get_max_threads
external omp_get_thread_num
integer omp_get_thread_num
external omp_get_num_procs
integer omp_get_num_procs
external omp_in_parallel
logical omp_in_parallel
external omp_set_dynamic
external omp_get_dynamic
logical omp_get_dynamic
external omp_get_cancellation
logical omp_get_cancellation
external omp_set_nested
external omp_get_nested
logical omp_get_nested
external omp_set_schedule
external omp_get_schedule
external omp_get_thread_limit
integer omp_get_thread_limit
external omp_set_max_active_levels
external omp_get_max_active_levels
integer omp_get_max_active_levels
external omp_get_level
integer omp_get_level
external omp_get_ancestor_thread_num
integer omp_get_ancestor_thread_num
external omp_get_team_size
integer omp_get_team_size
external omp_get_active_level
integer omp_get_active_level
external omp_set_default_device
external omp_get_default_device
integer omp_get_default_device
external omp_get_num_devices
integer omp_get_num_devices
external omp_get_num_teams
integer omp_get_num_teams
external omp_get_team_num
integer omp_get_team_num
external omp_is_initial_device
logical omp_is_initial_device
external omp_get_initial_device
integer omp_get_initial_device
external omp_get_max_task_priority
integer omp_get_max_task_priority

external omp_in_final
logical omp_in_final

integer ( omp_proc_bind_kind ) omp_get_proc_bind
external omp_get_proc_bind
integer omp_get_num_places
external omp_get_num_places
integer omp_get_place_num_procs
external omp_get_place_num_procs
external omp_get_place_proc_ids
integer omp_get_place_num
external omp_get_place_num
integer omp_get_partition_num_places
external omp_get_partition_num_places
external omp_get_partition_place_nums

external omp_init_lock
external omp_init_lock_with_hint
external omp_destroy_lock
external omp_set_lock
external omp_unset_lock
external omp_test_lock
logical omp_test_lock

external omp_init_nest_lock
external omp_init_nest_lock_with_hint
external omp_destroy_nest_lock
external omp_set_nest_lock
external omp_unset_nest_lock
external omp_test_nest_lock
integer omp_test_nest_lock

external omp_get_wtick
double precision omp_get_wtick
external omp_get_wtime
double precision omp_get_wtime
