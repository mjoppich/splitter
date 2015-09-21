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



namespace Ui {
class MainWindow;
}

class QDropEvent;
class QMimeData;

class ApplicationThread : public QThread
{
    Q_OBJECT

public:

    ~ApplicationThread()
    {

    }

protected:

    virtual void* performAction(void* pInData) = 0;

    void run()
    {

        // this thread does not need highest priority!
        QThread::currentThread()->setPriority(QThread::LowestPriority);

        m_pData = this->performAction(m_pInput);

        this->threadFinished(m_pData);
    }

    void* m_pInput;
    void* m_pData;

public slots:
    void threadFinished(void* pData)
    {
        emit getThreadFinished(pData);
    }
signals:
    void getThreadFinished(void* pData);

private:

};

class GTXthread : public ApplicationThread
{

public:
    GTXthread()
        : ApplicationThread()
    {

    }

    void* performAction(void* pInData)
    {

    }

};

class ExtendedBuffer : public QObject, public std::basic_streambuf<char> {
Q_OBJECT
public:

    ExtendedBuffer()
        : std::basic_streambuf<char>()
    {
    }

    ~ExtendedBuffer()
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
    void getApplicationFinished(void *pData);
    void receiveText(QString sString, QColor oColor);

    void on_pushButton_clicked();

protected:
    void dropEvent(QDropEvent *oDropEvent);
    void dragEnterEvent(QDragEnterEvent *oDragEnterEvent);

private:
    Ui::MainWindow *ui;

// controls

    bool m_bApplicationInProgress;
    QString m_sConfigFile;
    ExtendedBuffer* m_pBufferCout;
    ExtendedBuffer* m_pBufferCerr;
    ApplicationThread* m_pThread;
};

#endif // MAINWINDOW_H
