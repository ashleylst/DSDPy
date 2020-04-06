import threading
import time
from src import graph_processor as gp


class ProcessorThread(threading.Thread):

    def __init__(self, *args, **kwargs):
        super(ProcessorThread, self).__init__(*args, **kwargs)
        self.__flag = threading.Event()
        self.__flag.set()
        self.__running = threading.Event()
        self.__running.set()
        self.stopped = False
        self.lock = threading.Lock()

    def run(self):
        while not self._args[6][self._args[5]]:
            if self.stopped:
                return
            self.__flag.wait()
            self._args = gp.one_iteration(*self._args)
            print(time.time())
            time.sleep(1)
            if self.lock.locked():
                self.lock.release()
        print("graph processor finished.")

    def pause(self):
        print("is paused")
        self.__flag.clear()

    def resume(self):
        print("resume")
        self.__flag.set()

    def stop(self):
        print("stop")
        self.__flag.clear()
        self.__running.clear()
        self.stopped = True

    def get_lock(self):
        return self.lock

    def get_arg_info(self):
        return self._args

