import copy
import re
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)   # will display errors, but not warnings

class R_output(object):
    
    def capture_r_output(self, debug=False):
        """
        Capture and create a list of R console messages
        """

        # Record output #
        self.stdout = []
        self.stderr = []

        # Dummy functions #
        def add_to_stdout(line): self.stdout.append(line)
        def add_to_stderr(line): self.stderr.append(line)

        # Keep the old functions #
        self.stdout_orig = copy.deepcopy(rpy2.rinterface_lib.callbacks.consolewrite_print)
        self.stderr_orig = copy.deepcopy(rpy2.rinterface_lib.callbacks.consolewrite_warnerror)
        
        # If debug==False, uses python to print R lines after executing an R block 
        # If debug==True, will ugly print from R directly. Allows printing from R to troubleshoot errors.
        if not debug:

            # Set the call backs #
            rpy2.rinterface_lib.callbacks.consolewrite_print     = add_to_stdout
            rpy2.rinterface_lib.callbacks.consolewrite_warnerror = add_to_stderr

    
    def print_captured_r_output(self):
        """
        Cleans up R output for printing (e.g., removes "[1]" and end-of-line designators)
        """
        printable_lines = [line for line in self.stdout if line not in ['[1]', '\n']]
        printable_lines = [line for line in printable_lines if re.search("^\s*\[[0-9]+\]$", line) is None]
        printable_lines = [re.sub(r' \\n\"', "", line) for line in printable_lines]
        [print(line[2:-1]) for line in printable_lines]
        
        rpy2.rinterface_lib.callbacks.consolewrite_print     = self.stdout_orig
        rpy2.rinterface_lib.callbacks.consolewrite_warnerror = self.stderr_orig
