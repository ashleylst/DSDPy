import threading
import time
from PyQt5.QtCore import pyqtSignal, QThread
from src.basics import graph_processor as gp


class PProcessorThread(QThread):
    my_signal = pyqtSignal()

    def __init__(self, parent=None):
        super(QThread, self).__init__()
        self.stopped = False

    def run(self):
        self.my_signal.emit()
        pass


class ProcessorThread(threading.Thread):

    def __init__(self, event, *args, **kwargs):
        super(ProcessorThread, self).__init__(*args, **kwargs)
        self.pp = PProcessorThread()
        self.event = event
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
            # time.sleep(1)
            if self.lock.locked():
                self.lock.release()
        self.event.set()
        self.pp.start()
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

