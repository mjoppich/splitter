#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <iostream>
#include <fstream>
#include <omp.h>
#include <QMainWindow>
#include <QString>
#include <QTextEdit>
#include <QObject>
#include <QThread>

#include <libPTA/PTA_ErrorCorrect.h>

namespace Ui {
class MainWindow;
}

class QDropEvent;
class QMimeData;

class PTAErrorCorrectionThread : public QThread
{
    Q_OBJECT

public:
    void setConfigs(const char** pConfigs, unsigned int iConfigs)
    {
        m_iConfigs = iConfigs;
        m_pConfigs = (char**) pConfigs;
    }

    ~PTAErrorCorrectionThread()
    {
        free(m_pConfigs);
    }

protected:
    void run()
    {

        // this thread does not need highest priority!
        QThread::currentThread()->setPriority(QThread::LowestPriority);

        PTA::PTAStatistics* pStats = NULL;

        PTAErrorCorrect* pErrorCorrect = new PTAErrorCorrect(m_iConfigs, (const char**) m_pConfigs);
        pStats = pErrorCorrect->start();

        this->threadFinished(pStats);
    }

public slots:
    void threadFinished(PTA::PTAStatistics* pStats)
    {
        emit getPTAECFinished(pStats);
    }
signals:
    void getPTAECFinished(PTA::PTAStatistics* pStats);

private:
    int m_iConfigs;
    char** m_pConfigs;
};

class PTAExtendedBuffer : public QObject, public std::basic_streambuf<char> {
Q_OBJECT
public:

    PTAExtendedBuffer()
        : std::basic_streambuf<char>()
    {
    }

    ~PTAExtendedBuffer()
    {
    }

    void setTextColor(QColor oColor)
    {
        m_oColor = oColor;
    }

public slots:
    void transferText(QString sString)
    {
        emit sendText(sString, m_oColor);
    }
signals:
    void sendText(QString sString, QColor oColor);

protected:
    virtual int_type overflow(int_type v)
        {
            if (v == '\n')
            {
            }
            return v;
        }

    virtual std::streamsize xsputn(const char *p, std::streamsize n)
        {

            std::string sStdText(p, n);

            if ((p[n] == 4) || (p[n] == 0))
            {
                sStdText.append("\n");
            }

            QString sText(sStdText.c_str());

            this->transferText(sText);

            return n;
        }

private:
    QColor m_oColor;

};

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_oButtonStartEC_clicked();
    void on_oButtonConfig_clicked();
    void getPTAECFinished(PTA::PTAStatistics *pStats);
    void receiveText(QString sString, QColor oColor);

    void on_pushButton_clicked();

protected:
    void dropEvent(QDropEvent *oDropEvent);
    void dragEnterEvent(QDragEnterEvent *oDragEnterEvent);

private:
    Ui::MainWindow *ui;

// controls

    bool m_bCorrectionInProgress;
    QString m_sConfigFile;
    PTAExtendedBuffer* m_pBufferCout;
    PTAExtendedBuffer* m_pBufferCerr;
    PTAErrorCorrectionThread* m_pThread;
};

#endif // MAINWINDOW_H
