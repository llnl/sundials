#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------
# Send email notification if a SUNDIALS regression test status
# -----------------------------------------------------------------------------


def main():

    import argparse
    import os

    parser = argparse.ArgumentParser(
        description="Send email notification based on regression test status",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "teststatus",
        type=str,
        choices=["passed", "failed", "fixed"],
        help="Status of regression test",
    )

    parser.add_argument("testname", type=str, help="Name branch name or pull-request tested")

    parser.add_argument("testurl", type=str, help="URL for viewing test results")

    # parse command line args
    args = parser.parse_args()

    # name of regression test log file
    logfile = "suntest.log"

    # if log file exists add url, otherwise create log file
    if os.path.isfile(logfile):
        with open(logfile, "a") as log:
            log.write("View test output at:\n")
            log.write(args.testurl)
    else:
        print("Warning: log file not found")
        with open(logfile, "w") as log:
            log.write("Warning: log file not found\n")
            log.write("View test output at:\n")
            log.write(args.testurl)

    # determine notification recipient
    special_branches = ["main", "develop", "release"]

    if any(branch in args.testname for branch in special_branches):
        # SUNDIALS developers list
        recipient = "sundials-devs@llnl.gov"
    else:
        # author of most recent commit
        cmd = "git log --format='%ae' -1"
        recipient = runCommand(cmd).rstrip().decode("UTF-8")

        # check if the last commit was a CI merge
        if recipient == "nobody@nowhere":
            cmd = "git log HEAD~1 --pretty=format:'%ae' -1"
            recipient = runCommand(cmd).rstrip().decode("UTF-8")

    # send notification if tests fail, log file not found, or fixed
    if args.teststatus == "failed":

        subject = "FAILED: SUNDIALS " + args.testname + " failed regression tests"
        print("Tests failed, sending notification to", recipient)
        sendEmail(recipient, subject, logfile)

    elif args.teststatus == "fixed":

        subject = "FIXED: SUNDIALS " + args.testname + " passed regression tests"
        print("Tests fixed, sending notification to", recipient)
        sendEmail(recipient, subject, logfile)

    else:
        print("Tests passed, no notifications sent")


#
# run external command
#
def runCommand(cmd):

    import subprocess

    cmdout = subprocess.check_output(cmd, shell=True)

    return cmdout


#
# send email
#
def sendEmail(recipient, subject, message):

    import smtplib
    from email.message import EmailMessage

    # Open a plain text file for reading. Assumed that
    # the text file contains only ASCII characters.
    with open(message) as fp:
        # Create a text/plain message
        msg = EmailMessage()
        msg.set_content(fp.read())

    # sender's email address
    sender = "SUNDIALS.suntest@llnl.gov"

    # email settings
    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = recipient

    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP("smtp.llnl.gov")
    s.send_message(msg)
    s.quit()


#
# just run the main routine
#
if __name__ == "__main__":
    main()
