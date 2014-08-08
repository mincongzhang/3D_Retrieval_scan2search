// MeshRetrievalMFCDlg.h : header file
//

#pragma once

#include "OpenGLControl.h"
#include "MeshOperation.h"

// CMeshRetrievalMFCDlg dialog
class CMeshRetrievalMFCDlg : public CDialog
{
	// Construction
public:
	CMeshRetrievalMFCDlg(CWnd* pParent = NULL);	// standard constructor

	COpenGLControl m_oglWindow;

	// Dialog Data
	enum { IDD = IDD_MeshRetrievalMFC_DIALOG };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


	// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();

	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnTimer(UINT_PTR nIDEvent);
	afx_msg void OnBnClickedLoad();
	afx_msg void OnBnClickedNoise();
	afx_msg void OnBnClickedNormalize();
	afx_msg void OnBnClickedSketch();
	afx_msg void OnBnClickedSpharm();
	afx_msg void OnBnClickedCandidate1();
	afx_msg void OnBnClickedCandidate2();
	afx_msg void OnBnClickedCandidate3();
	afx_msg void OnBnClickedCandidate4();
	afx_msg void OnBnClickedCandidate5();
	afx_msg void OnBnClickedCandidate6();
	afx_msg void OnBnClickedRasterize();
	afx_msg void OnBnClickedBatchtransform();
	afx_msg void OnBnClickedDenoise();
	afx_msg void OnBnClickedRetrieve();
	afx_msg void OnBnClickedX();
	afx_msg void OnBnClickedY();
	afx_msg void OnBnClickedZ();
};
