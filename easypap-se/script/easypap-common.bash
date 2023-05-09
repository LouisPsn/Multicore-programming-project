#!/usr/bin/env bash

EASYPAPDIR=${EASYPAPDIR:-.}

TRACEDIR=${TRACEDIR:-${EASYPAPDIR}/traces/data}
TRACEFILE=${TRACEDIR}/ezv_trace_current.evt
OLDTRACEFILE=${TRACEDIR}/ezv_trace_previous.evt

SIMU=${EASYPAPDIR}/bin/easypap
VIEW=${EASYPAPDIR}/traces/bin/easyview
