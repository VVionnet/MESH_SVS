
integer omp_lock_kind
integer omp_nest_lock_kind
integer omp_lock_hint_kind
! this selects an integer that is large enough to hold a 64 bit integer
parameter ( omp_lock_kind = selected_int_kind( 10 ) )
parameter ( omp_nest_lock_kind = selected_int_kind( 10 ) )
parameter ( omp_lock_hint_kind = selected_int_kind( 10 ) )

integer omp_sched_kind
! this selects an integer that is large enough to hold a 32 bit integer
parameter ( omp_sched_kind = selected_int_kind( 8 ) )
integer ( omp_sched_kind ) omp_sched_static
parameter ( omp_sched_static = 1 )
integer ( omp_sched_kind ) omp_sched_dynamic
parameter ( omp_sched_dynamic = 2 )
integer ( omp_sched_kind ) omp_sched_guided
parameter ( omp_sched_guided = 3 )
integer ( omp_sched_kind ) omp_sched_auto
parameter ( omp_sched_auto = 4 )

integer omp_proc_bind_kind
parameter ( omp_proc_bind_kind = selected_int_kind( 8 ) )
integer ( omp_proc_bind_kind ) omp_proc_bind_false
parameter ( omp_proc_bind_false = 0 )
integer ( omp_proc_bind_kind ) omp_proc_bind_true
parameter ( omp_proc_bind_true = 1 )
integer ( omp_proc_bind_kind ) omp_proc_bind_master
parameter ( omp_proc_bind_master = 2 )
integer ( omp_proc_bind_kind ) omp_proc_bind_close
parameter ( omp_proc_bind_close = 3 )
integer ( omp_proc_bind_kind ) omp_proc_bind_spread
parameter ( omp_proc_bind_spread = 4 )

integer ( omp_lock_hint_kind ) omp_lock_hint_none
parameter ( omp_lock_hint_none = 0 )
integer ( omp_lock_hint_kind ) omp_lock_hint_uncontended
parameter ( omp_lock_hint_uncontended = 1 )
integer ( omp_lock_hint_kind ) omp_lock_hint_contended
parameter ( omp_lock_hint_contended = 2 )
integer ( omp_lock_hint_kind ) omp_lock_hint_nonspeculative
parameter ( omp_lock_hint_nonspeculative = 4 )
integer ( omp_lock_hint_kind ) omp_lock_hint_speculative
parameter ( omp_lock_hint_speculative = 8 )
