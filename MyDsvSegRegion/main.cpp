#define _CRT_SECURE_NO_WARNINGS
#include "define.h"

TRANSINFO	calibInfo;

HANDLE	dfp = INVALID_HANDLE_VALUE;
int		dsbytesiz = sizeof(point3d) * 2 + sizeof(ONEVDNDATA);
int		dFrmNum = 0;
int		dFrmNo = 0;

RMAP	rm;
DMAP	dm;
DMAP	gm, ggm;

ONEDSVFRAME	*onefrm;

bool LoadCalibFile(char *szFile)
{
	char			i_line[200];
	FILE			*fp;
	MATRIX			rt;

	fp = fopen(szFile, "r");
	if (!fp)
		return false;

	rMatrixInit(calibInfo.rot);

	int	i = 0;
	while (1) {
		if (fgets(i_line, 80, fp) == NULL)
			break;

		if (_strnicmp(i_line, "rot", 3) == 0) {
			strtok(i_line, " ,\t\n");
			calibInfo.ang.x = atof(strtok(NULL, " ,\t\n"))*topi;
			calibInfo.ang.y = atof(strtok(NULL, " ,\t\n"))*topi;
			calibInfo.ang.z = atof(strtok(NULL, " ,\t\n"))*topi;
			createRotMatrix_ZYX(rt, calibInfo.ang.x, calibInfo.ang.y, calibInfo.ang.z);
			rMatrixmulti(calibInfo.rot, rt);
			continue;
		}

		if (_strnicmp(i_line, "shv", 3) == 0) {
			strtok(i_line, " ,\t\n");
			calibInfo.shv.x = atof(strtok(NULL, " ,\t\n"));
			calibInfo.shv.y = atof(strtok(NULL, " ,\t\n"));
			calibInfo.shv.z = atof(strtok(NULL, " ,\t\n"));
		}
	}
	fclose(fp);

	return true;
}

void SmoothingData()
{
	int maxcnt = 3;

	for (int y = 0; y<rm.len; y++) {
		for (int x = 1; x<(rm.wid - 1); x++) {
			if (rm.pts[y*rm.wid + (x - 1)].i && !rm.pts[y*rm.wid + x].i) {

				int xx;
				for (xx = x + 1; xx<rm.wid; xx++) {
					if (rm.pts[y*rm.wid + xx].i)
						break;
				}
				if (xx >= rm.wid)
					continue;
				int cnt = xx - x + 1;
				if (cnt>maxcnt) {
					x = xx;
					continue;
				}
				point3fi *p1 = &rm.pts[y*rm.wid + (x - 1)];
				point3fi *p2 = &rm.pts[y*rm.wid + xx];
				double dis = ppDistance3fi(p1, p2);
				double rng = max(p2r(p1), p2r(p2));
				double maxdis = min(MAXSMOOTHERR, max(BASEERROR, HORIERRFACTOR*cnt*rng));
				if (dis<maxdis) {
					for (int xxx = x; xxx<xx; xxx++) {
						point3fi *p = &rm.pts[y*rm.wid + xxx];
						p->x = (p2->x - p1->x) / cnt * (xxx - x + 1) + p1->x;
						p->y = (p2->y - p1->y) / cnt * (xxx - x + 1) + p1->y;
						p->z = (p2->z - p1->z) / cnt * (xxx - x + 1) + p1->z;
						p->i = 1;
					}
				}
				x = xx;
			}
		}
	}
}

void CorrectPoints()
{
	MAT2D	rot1, rot2;

	//transform points to the vehicle frame of onefrm->dsv[0]
	//src: block i; tar: block 0

	//rot2: R_tar^{-1}
	rot2[0][0] = cos(-onefrm->dsv[0].ang.z);
	rot2[0][1] = -sin(-onefrm->dsv[0].ang.z);
	rot2[1][0] = sin(-onefrm->dsv[0].ang.z);
	rot2[1][1] = cos(-onefrm->dsv[0].ang.z);

	for (int i = 1; i<BKNUM_PER_FRM; i++) {
		for (int j = 0; j<PTNUM_PER_BLK; j++) {
			if (!onefrm->dsv[i].points[j].i)
				continue;

			rotatePoint3fi(onefrm->dsv[i].points[j], calibInfo.rot);
			shiftPoint3fi(onefrm->dsv[i].points[j], calibInfo.shv);
			rotatePoint3fi(onefrm->dsv[i].points[j], onefrm->dsv[i].rot);

			//rot1: R_tar^{-1}*R_src
			rot1[0][0] = cos(onefrm->dsv[i].ang.z - onefrm->dsv[0].ang.z);
			rot1[0][1] = -sin(onefrm->dsv[i].ang.z - onefrm->dsv[0].ang.z);
			rot1[1][0] = sin(onefrm->dsv[i].ang.z - onefrm->dsv[0].ang.z);
			rot1[1][1] = cos(onefrm->dsv[i].ang.z - onefrm->dsv[0].ang.z);


			//shv: SHV_src-SHV_tar
			point2d shv;
			shv.x = onefrm->dsv[i].shv.x - onefrm->dsv[0].shv.x;
			shv.y = onefrm->dsv[i].shv.y - onefrm->dsv[0].shv.y;

			point2d pp;
			pp.x = onefrm->dsv[i].points[j].x; pp.y = onefrm->dsv[i].points[j].y;
			rotatePoint2d(pp, rot1);	//R_tar^{-1}*R_src*p
			rotatePoint2d(shv, rot2);	//R_tar^{-1}*(SHV_src-SHV_tar)
			shiftPoint2d(pp, shv);		//p'=R_tar^{-1}*R_src*p+R_tar^{-1}*(SHV_src-SHV_tar)
			onefrm->dsv[i].points[j].x = pp.x;
			onefrm->dsv[i].points[j].y = pp.y;
		}
	}

	for (int ry = 0; ry<rm.len; ry++) {
		for (int rx = 0; rx<rm.wid; rx++) {
			int i = rm.idx[ry*rm.wid + rx].x;
			int j = rm.idx[ry*rm.wid + rx].y;
			if (!i && !j)
				continue;
			rm.pts[ry*rm.wid + rx] = onefrm->dsv[i].points[j];
		}
	}
}

void ProcessOneFrame()
{
	//���ɾ���ͼ��֡
	GenerateRangeView();

	//����calib�����������ת������������ϵ�����ݳ���Ƕ�roll��pitch��������֡��ˮƽ�����ݵ�ת������0�����ݰ�����Ǻ�λ������Ӧ�ĳ�������ϵ
	CorrectPoints();

	//��ÿһ�������ж�����Ч����㣨cnt<5��Լˮƽ1��)�����ڲ岹�룬������Щ��Ч�㴦�ᱻ��Ϊ�Ǳ߽��
	SmoothingData();

	//�ָ�
	//��һ������ע�߽��ContourExtraction();
	//�ڶ���������������ʽ��ע�����ڵ�RegionGrow()
	memset(rm.regionID, 0, sizeof(int)*rm.wid*rm.len);
	rm.regnum = 0;
	ContourSegger();

	//Ϊÿ����������һ��segbuf�����ڷ��ࡢĿǰ����ȡ����������
	if (rm.regnum) {
		rm.segbuf = new SEGBUF[rm.regnum];
		memset(rm.segbuf, 0, sizeof(SEGBUF)*rm.regnum);
		Region2Seg();
	}

	//���ɿ��ӻ�����ͼ������
	//DrawRangeView ();

	//��ȫ��DEMת������ǰ��������ϵ��
	//PredictGloDem(gm, ggm);
	//printf("PredictGloDem completed\n");
	//���ɵ�֡���ݵ�DEM
	GenerateLocDem(dm);
	//printf("GenerateLocDem completed\n");

	//�õ�ǰ֡DEM����ȫ��DEM
	//UpdateGloDem(gm, dm);
	
	//��ȡ��·������
	//ExtractRoadCenterline(gm);

	//������温���ͺ���ǣ��������
	//LabelRoadSurface(gm);

	//��ȡ·���ϵ��ϰ���
	//LabelObstacle(gm);

	//���ɿ��ӻ���֡����DEM
	//DrawDem(dm);

	//���ɿ��ӻ�ȫ��DEM
	//DrawDem(gm);
	
	if (rm.segbuf)
		delete[]rm.segbuf;

}

//��ȡһ֡vel32���ݣ�һ֡Ϊ180��12��32������㣩���浽onefrm->dsv��δ������ת��
BOOL ReadOneDsvFrame()
{
	DWORD	dwReadBytes;
	int		i;

	for (i = 0; i<BKNUM_PER_FRM; i++) {
		if (!ReadFile(dfp, (ONEDSVDATA *)&onefrm->dsv[i], dsbytesiz, &dwReadBytes, NULL) || (dsbytesiz != dwReadBytes))
			break;
		//		createRotMatrix_ZYX(onefrm->dsv[i].rot, onefrm->dsv[i].ang.x, onefrm->dsv[i].ang.y , onefrm->dsv[i].ang.z ) ; 
		createRotMatrix_ZYX(onefrm->dsv[i].rot, onefrm->dsv[i].ang.x, onefrm->dsv[i].ang.y, 0);
	}

	if (i<BKNUM_PER_FRM)
		return FALSE;
	else
		return TRUE;
}

void CallbackLocDem(int event, int x, int y, int flags, void *ustc)
{
	static CvPoint lu, rb;

	if (event == CV_EVENT_LBUTTONDOWN) {
		lu.x = x; lu.y = y;
	}
	if (event == CV_EVENT_LBUTTONUP) {

		rb.x = x; rb.y = y;
		IplImage *tmp = cvCreateImage(cvSize(dm.wid, dm.len), IPL_DEPTH_8U, 3);
		cvCopy(dm.lmap, tmp);
		cvRectangle(dm.lmap, lu, rb, cvScalar(255, 255, 0), 3);
		cvShowImage("ldemlab", dm.lmap);

		int ix, iy;
		for (iy = min(lu.y, rb.y); iy <= max(lu.y, rb.y); iy++)
			for (ix = min(lu.x, rb.x); ix <= max(lu.x, rb.x); ix++)
				printf("%d, %d, %.3f,%.3f\n", ix, iy, dm.demg[iy*dm.wid + ix], dm.demhmin[iy*dm.wid + ix]);
		cvReleaseImage(&tmp);
	}
}

void JumpTo(int fno)
{
	LARGE_INTEGER li;

	li.QuadPart = __int64(dsbytesiz)*BKNUM_PER_FRM*fno;
	li.LowPart = SetFilePointer(dfp, li.LowPart, &li.HighPart, FILE_BEGIN);
	dFrmNo = fno;
}

//���������
void DoProcessing()
{

	DWORD dwSizeLow, dwSizeHigh, dwError;  //��λ����λ���������
	dwSizeLow = GetFileSize((HANDLE)dfp, &dwSizeHigh);  //��������
	if (dwSizeLow == 0xFFFFFFFF && (dwError = GetLastError()) != NO_ERROR) {
		return;
	}
	else {
		LONGLONG llSize, llPow;
		llPow = 4294967296; //(LONGLONG)pow(2,32);
		llSize = dwSizeHigh * llPow + dwSizeLow; //�����Ǵ���4G�Ĺ���
		dFrmNum = llSize / 180 / dsbytesiz;
		printf("Total Frame Number is %d\n", dFrmNum);
	}

	SetFilePointer(dfp, 0, 0, FILE_BEGIN);

	InitRmap(&rm);
	InitDmap(&dm);
	//InitDmap(&gm);
	//InitDmap(&ggm);
	onefrm = new ONEDSVFRAME[1];

	//IplImage * col = cvCreateImage(cvSize(1024, rm.len * 3), IPL_DEPTH_8U, 3);

	//CvFont font;
	//cvInitFont(&font, CV_FONT_HERSHEY_DUPLEX, 1, 1, 0, 2);

	int waitkeydelay = 0;

	//	JumpTo (dFrmNum/2);
	dFrmNo = 0;

	while (ReadOneDsvFrame())
	{
		if (dFrmNo % 100 == 0)
			printf("%d (%d)\n", dFrmNo, dFrmNum);

		//ÿһ֡�Ĵ���
		ProcessOneFrame();

		SaveDEM(dm, dFrmNo, dFrmNum);
		//���ӻ�
		//		cvResize (rm.rMap, col);
		/*
		char str[10];
		sprintf(str, "Fno%d", dFrmNo);
		cvPutText(dm.lmap, str, cvPoint(50, 50), &font, CV_RGB(0, 0, 255));
		*/
		//		cvShowImage("range image",col);
		//		cvResize (rm.lMap, col);
		//		cvShowImage("region",col);
		/*
		if (dm.lmap) cvShowImage("ldemlab", dm.lmap);
		if (gm.lmap) cvShowImage("gdemlab", gm.lmap);
		if (gm.smap) cvShowImage("gsublab", gm.smap);

		cv::setMouseCallback("gsublab", CallbackLocDem, 0);

		char WaitKey;
		WaitKey = cvWaitKey(waitkeydelay);
		if (WaitKey == 27)
			break;
		if (WaitKey == 'z') {
			if (waitkeydelay == 1)
				waitkeydelay = 0;
			else
				waitkeydelay = 1;
		}
		*/
		dFrmNo++;
	}

	ReleaseRmap(&rm);
	ReleaseDmap(&dm);
	//ReleaseDmap(&gm);
	//ReleaseDmap(&ggm);
	//cvReleaseImage(&col);
	delete[]onefrm;
}

void main(int argc, char *argv[])
{

	if (argc<3) {
		fprintf(stderr, "Usage : %s [infile] [calibfile]\n", argv[0]);
		fprintf(stderr, "[infile] DSV file.\n");
		fprintf(stderr, "[calibfile] define the calibration parameters of the DSV file.\n");
		fprintf(stderr, "[outfile] segmentation results to DSVL file.\n");
		fprintf(stderr, "[seglog] data association results to LOG file.\n");
		fprintf(stderr, "[videooutflg] 1: output video to default files, 0: no output.\n");
		exit(1);
	}

	if (!LoadCalibFile(argv[2])) {
		fprintf(stderr, "Invalid calibration file : %s.\n", argv[2]);
		getchar();
		exit(1);
	}

	dfp = CreateFile(argv[1], GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
	if (dfp == INVALID_HANDLE_VALUE) {
		printf("File open failure : %s\n", argv[1]);
		getchar();
		exit(1);
	}

	DoProcessing();

	fprintf(stderr, "Labeling succeeded.\n");

	CloseHandle(dfp);

	system("pause");
	exit(1);
}

