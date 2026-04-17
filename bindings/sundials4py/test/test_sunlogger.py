#!/usr/bin/env python3

import pytest
import sys
from io import StringIO
from sundials4py.core import *


def test_custom_queue_and_flush_and_restore_default():
    """
    Test that custom queue/flush functions work and that passing
    None values restores the default behavior (errors to stderr).
    """
    # Create logger
    err, logger = SUNLogger_Create(SUN_COMM_NULL, 0)
    assert err == SUN_SUCCESS
    assert logger is not None

    # Storage for custom handler messages
    messages = []

    # Custom queue function that stores messages in a list
    def custom_queue_fn(logger, lvl, prefix, rank, scope, label, payload, content):
        msg = f"[{prefix}][rank {rank}][{scope}][{label}] {payload}"
        messages.append(msg)
        return SUN_SUCCESS

    # Custom flush function
    def custom_flush_fn(logger, lvl, content):
        return SUN_SUCCESS

    # Set custom queue and flush functions
    status = SUNLogger_SetQueueAndFlushMsgFns(logger, custom_queue_fn, custom_flush_fn)
    assert status == SUN_SUCCESS

    # Queue and flush a message with custom handlers
    status = SUNLogger_QueueMsg(
        logger, SUN_LOGLEVEL_ERROR, "test_scope", "test_label", "custom handler message"
    )
    assert status == SUN_SUCCESS

    status = SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR)
    assert status == SUN_SUCCESS

    # Verify the custom handler was called
    assert len(messages) == 1
    assert "custom handler message" in messages[0]
    assert "test_scope" in messages[0]
    assert "test_label" in messages[0]

    # Now restore default behavior by passing None values
    status = SUNLogger_SetQueueAndFlushMsgFns(logger, None, None)
    assert status == SUN_SUCCESS

    # Capture stderr to verify default behavior
    old_stderr = sys.stderr
    sys.stderr = StringIO()

    # Queue and flush another error message
    status = SUNLogger_QueueMsg(
        logger, SUN_LOGLEVEL_ERROR, "default_scope", "default_label", "default error message"
    )
    assert status == SUN_SUCCESS

    status = SUNLogger_Flush(logger, SUN_LOGLEVEL_ERROR)
    assert status == SUN_SUCCESS

    # Get stderr output and restore original stderr
    stderr_output = sys.stderr.getvalue()
    sys.stderr = old_stderr

    # Verify the message appears on stderr (default behavior)
    assert stderr_output.count("\n") == 1
    assert "default error message" in stderr_output

    # Verify the custom handler was NOT called again
    assert len(messages) == 1
