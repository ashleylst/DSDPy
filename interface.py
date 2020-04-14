import sys
import os
import threading
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
import matplotlib.pyplot as plt
from src.util import processorthread, model, uiwindow
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from src.basics import output, graph_processor


class MyWindow(QMainWindow, uiwindow.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MyWindow, self).__init__(parent)
        self.model = model.Model()
        self.procthread = None
        self.lock = None
        self.simuarg = None
        self.simumode = 'bng'
        self.paused = False
        self.canvas = None
        self.setupUi(self)

    def refresh(self):
        self.lineEdit.setText(self.model.get_filename())
        self.textEdit.setText(self.model.get_filecontents())

    def returnPressedSlot(self):
        fname = self.lineEdit.text()
        if self.model.is_valid(fname):
            self.model.set_filename(self.lineEdit.text())
            self.refresh()
        else:
            m = QtWidgets.QMessageBox()
            m.setText("Invalid file name!\n" + fname)
            m.setIcon(QtWidgets.QMessageBox.Warning)
            m.setStandardButtons(QtWidgets.QMessageBox.Ok
                                 | QtWidgets.QMessageBox.Cancel)
            m.setDefaultButton(QtWidgets.QMessageBox.Cancel)
            ret = m.exec_()
            self.lineEdit.setText("")
            self.refresh()

    def browseSlot(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fname, filetype = QFileDialog.getOpenFileName(self,
                                                          "Select File",
                                                          "./",
                                                          "All Files (*);;Text Files (*.txt)")
        if fname:
            self.debugPrint("Reading File.")
            self.model.set_filename(fname)
            self.refresh()

    def runSlot(self):
        pass

    def in_generation(self):
        self.pushButton_pause.setEnabled(True)
        self.pushButton_stop.setEnabled(True)

        self.pushButton_analyze.setEnabled(False)
        self.pushButton_show.setEnabled(False)
        self.pushButton_simu.setEnabled(False)
        self.comboBox_simumode.setEnabled(False)

    def out_of_generation(self):
        self.pushButton_pause.setEnabled(False)
        self.pushButton_stop.setEnabled(False)

        self.pushButton_analyze.setEnabled(True)
        self.pushButton_show.setEnabled(True)
        self.pushButton_simu.setEnabled(True)
        self.comboBox_simumode.setEnabled(True)

    def my_event(self):
        self.debugPrint("Generation Finished.")
        self.out_of_generation()
        self.pushButton_simu.setEnabled(True)
        self.comboBox_simumode.setEnabled(True)
        self.showSlot()

    def submitSlot(self):
        if not self.procthread:
            return
        specieslist, speciesidmap, reactionlist, kinetics, indexlist, cursor, visited = self.procthread.get_arg_info()
        initnames, concentrations, outdir, simupara, initlen = self.simuarg
        try:
            x, y, obs = graph_processor.simulation(specieslist, reactionlist, initlen, initnames,
                                                               concentrations, outdir, simupara, self.simumode)
        except:
            self.debugPrint('Simulation Error.')
            return
        self.display_output_img(x, y, obs, option=self.simumode)

    def analyzeSlot(self):
        # Take care of submitting more than once
        if self.procthread:
            if self.procthread.is_alive():
                self.procthread.stop()

        text = self.textEdit.toPlainText()

        if self.model.get_filename() is None:
            if text == '':
                m = QtWidgets.QMessageBox()
                m.setText("Please provide an input to submit.\n")
                m.setIcon(QtWidgets.QMessageBox.Warning)
                m.setStandardButtons(QtWidgets.QMessageBox.Ok)
                m.setDefaultButton(QtWidgets.QMessageBox.Ok)
                ret = m.exec_()
                return
            else:
                m = QtWidgets.QMessageBox()
                m.setText("Do you want to save the input?\n")
                m.setIcon(QtWidgets.QMessageBox.Warning)
                m.setStandardButtons(QtWidgets.QMessageBox.Ok
                                     | QtWidgets.QMessageBox.Cancel)
                m.setDefaultButton(QtWidgets.QMessageBox.Ok)
                ret = m.exec_()
                if ret == m.Ok:
                    try:
                        os.mkdir('./input')
                    except FileExistsError:
                        pass
                    file = open('./input/your_input.txt', 'w+')
                    file.write(text)
                    file.close()
                    self.model.set_filename('./input/your_input.txt')

        elif self.model.get_filecontents() != text:
            m = QtWidgets.QMessageBox()
            m.setText("Do you want to save the change to input file?\n")
            m.setIcon(QtWidgets.QMessageBox.Warning)
            m.setStandardButtons(QtWidgets.QMessageBox.Ok
                                 | QtWidgets.QMessageBox.Cancel)
            m.setDefaultButton(QtWidgets.QMessageBox.Ok)
            ret = m.exec_()
            if ret == m.Ok:
                file = open(self.model.fileName, 'w+')
                file.write(self.textEdit.toPlainText())
                file.close()
                self.model.set_filename(self.model.fileName)

        self.debugPrint("Submitted File.")
        try:
            info, initnames, concentrations, outdir, simupara, initlen = graph_processor.initiation(text=text)
        except:
            self.debugPrint('Input Error.')
            return
        self.simuarg = (initnames, concentrations, outdir, simupara, initlen)

        event = threading.Event()
        self.procthread = processorthread.ProcessorThread(event, args=info)
        self.procthread.setDaemon(True)
        self.lock = self.procthread.get_lock()
        self.paused = False
        self.in_generation()

        self.procthread.start()

        # show results when finished
        self.procthread.pp.my_signal.connect(self.my_event)

    def pauseSlot(self):
        if not self.procthread:
            return
        if self.paused:
            self.resumeSlot()
            return

        print("need to pause")
        self.lock.acquire()
        self.procthread.pause()
        self.pushButton_pause.setText('Resume')
        self.pushButton_show.setEnabled(True)
        self.paused = True

    def resumeSlot(self):
        if not self.procthread:
            return
        print("need to resume")
        if self.lock.locked():
            self.lock.release()
        self.procthread.resume()
        self.paused = False
        self.pushButton_pause.setText("Pause")
        self.pushButton_show.setEnabled(False)

    def stopSlot(self):
        if not self.procthread:
            return
        if not self.paused:
            self.lock.acquire()
        self.pushButton_pause.setText("Pause")
        self.out_of_generation()
        self.procthread.stop()

    def showSlot(self):
        if not self.procthread:
            return
        specieslist, speciesidmap, reactionlist, kinetics, indexlist, cursor, visited = self.procthread.get_arg_info()
        try:
            text = graph_processor.post_enumeration(specieslist, reactionlist)
        except:
            self.debugPrint('Generation Error.')
            return
        self.display_output_txt(text)
        self.pushButton_save.setEnabled(True)

    def saveSlot(self):
        text = self.debugTextBrowser.toPlainText()
        if len(text) == 0:
            return

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fname, _ = QFileDialog.getSaveFileName(self, "QFileDialog.getSaveFileName()", "",
                                                  "All Files (*);;Text Files (*.txt)", options=options)
        if fname:
            file = open(fname, 'w+')
            file.write(text)
            file.close()

    def simuSlot(self):
        mode = self.comboBox_simumode.currentText()
        if mode == 'Stochastic':
            self.simumode = 'bng'
        else:
            self.simumode = 'scipy'
        self.debugPrint('Simulation mode set to ' + mode)

    def debugPrint(self, msg):
        self.debugTextBrowser.append(msg)

    def display_output_txt(self, text):
        self.reneTextBrowser.setText(text)
        # self.textBrowser.append(open(fname, 'r').read())

    def display_output_img(self, x, y, obs, colormap='tab10', option='bng'):
        obslen = len(obs)
        if self.canvas is None:
            sc = output.Canvas(self, width=5, height=6, dpi=100)
            self.canvas = sc
            toolbar = NavigationToolbar(sc, self)

            layout = QtWidgets.QVBoxLayout()
            layout.addWidget(toolbar)
            layout.addWidget(sc)

            # Create a placeholder widget to hold our toolbar and canvas.
            widget = QtWidgets.QWidget(self.frame_4)
            widget.setLayout(layout)
            self.verticalLayout_3.addWidget(widget)
        else:
            sc = self.canvas
            sc.axes.cla()

        cmap = plt.cm.get_cmap(colormap, obslen)

        for i in range(0, obslen):
            label = obs[i].name[3:]
            if option == 'bng':
                sc.axes.plot(x, y[:, i], label=label, c=cmap(i))
            elif option == 'scipy':
                sc.axes.plot(x, y[obs[i].name], label=label, c=cmap(i))

        sc.axes.set_xlabel("Time (s)")
        sc.axes.set_ylabel("Complexes")
        sc.axes.legend(bbox_to_anchor=(1.01, 1))
        sc.axes.set_xmargin(0.5)

        sc.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    myWin = MyWindow()
    myWin.setWindowTitle('DSDPy')
    myWin.show()
    sys.exit(app.exec_())
