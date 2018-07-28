﻿#include "define.h"


void DrawRangeView()
{
	int x, y;

	cvZero(rm.rMap);	//����ͼ��
	cvZero(rm.lMap);	//�ָ�ͼ��
						//���ָ������������ӻ�λͼ��ɫ
	for (y = 0; y<rm.len; y++) {
		for (x = 0; x<rm.wid; x++) {
			if (!rm.pts[y*rm.wid + x].i) {
				rm.rMap->imageData[(y*rm.wid + x) * 3] = 255;
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 1] = 255;
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 2] = 255;
				continue;
			}
			if (rm.pts[y*rm.wid + x].z<0.0) {
				rm.rMap->imageData[(y*rm.wid + x) * 3] = 0;
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 1] = 0;
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 2] = BOUND(nint(-rm.pts[y*rm.wid + x].z*100.0), 0, 255);
			}
			else {
				rm.rMap->imageData[(y*rm.wid + x) * 3] = BOUND(nint(rm.pts[y*rm.wid + x].z*100.0), 0, 255);
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 1] = 0;
				rm.rMap->imageData[(y*rm.wid + x) * 3 + 2] = 0;
			}
			if (rm.regionID[y*rm.wid + x] == EDGEPT) {
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 2] = 128;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 1] = 128;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 0] = 255;
			}
			else if (rm.regionID[y*rm.wid + x] == NONVALID)
			{
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 2] = 255;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 1] = 255;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 0] = 255;
			}
			else if (rm.regionID[y*rm.wid + x] == UNKNOWN)
			{
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 2] = 64;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 1] = 64;
				rm.lMap->imageData[(y*rm.wid + x) * 3 + 0] = 64;
			}
			else {
				SEGBUF *segbuf = &rm.segbuf[rm.regionID[y*rm.wid + x]];
				if (segbuf->ptnum) {
					rm.lMap->imageData[(y*rm.wid + x) * 3 + 0] = nint(fabs(segbuf->norm.x)*255.0);
					rm.lMap->imageData[(y*rm.wid + x) * 3 + 1] = nint(fabs(segbuf->norm.y)*255.0);
					rm.lMap->imageData[(y*rm.wid + x) * 3 + 2] = nint(segbuf->norm.z*255.0);
				}
			}
		}
	}
}

void GenerateRangeView()
{
	memset(rm.pts, 0, sizeof(point3fi)*rm.wid*rm.len);	//����ͼ���Ӧ�ļ��������
	memset(rm.idx, 0, sizeof(point2i)*rm.wid*rm.len);
	memset(rm.di, 0, sizeof(BYTE)*rm.wid*rm.len);

	//���ɾ���ͼ����Ӧ������
	for (int i = 0; i<BKNUM_PER_FRM; i++) {
		for (int j = 0; j<LINES_PER_BLK; j++) {
			for (int k = 0; k<PNTS_PER_LINE; k++) {
				point3fi *p = &onefrm->dsv[i].points[j*PNTS_PER_LINE + k];
				if (!p->i)
					continue;

				float rng = sqrt(sqr(p->x) + sqr(p->y) + sqr(p->z));
				float angv = asin(p->z / rng);
				float angh = atan2(p->y, p->x);
				float ix = (angh - rm.h0) / rm.hres;
				float iy = (angv - rm.v0) / rm.vres;

				int x0, y0, x1, y1;
				int inten;
				x0 = int(ix); y0 = int(iy);
				x1 = int(ix) + 1; y1 = int(iy) + 1;
				for (int y = y0; y <= y1; y++) {
					if (y<0 || y >= rm.len) continue;
					for (int x = x0; x <= x1; x++) {
						if (x<0 || x >= rm.wid) continue;
						if (y == 0 && x>416)
							y = y;
						inten = sqrt(sqr(x - ix) + sqr(y - iy)) * 100 + 10;
						if (!rm.pts[y*rm.wid + x].i) {
							//							rm.pts[y*rm.wid+x] = onefrm->dsv[i].points[j*PNTS_PER_LINE+k];
							rm.idx[y*rm.wid + x].x = i;
							rm.idx[y*rm.wid + x].y = j * PNTS_PER_LINE + k;
							rm.di[y*rm.wid + x] = inten;
						}
						else {
							if (inten>rm.di[y*rm.wid + x]) {
								//								rm.pts[y*rm.wid+x] = onefrm->dsv[i].points[j*PNTS_PER_LINE+k];
								rm.idx[y*rm.wid + x].x = i;
								rm.idx[y*rm.wid + x].y = j * PNTS_PER_LINE + k;
								rm.di[y*rm.wid + x] = inten;
							}
						}
					}
				}
			}
		}
	}
}

void InitRmap(RMAP *rm)
{
	rm->wid = LINES_PER_BLK * BKNUM_PER_FRM / 2;
	rm->wid /= 2.0;								//downsample horizontal scans to 1/2
	rm->len = PNTS_PER_LINE * 2;
	rm->hres = M_PI * 2.0 / (rm->wid - 1);
	rm->vres = (VMAXANG - VMINANG) / (rm->len - 1);
	rm->h0 = -M_PI;
	rm->v0 = VMINANG;
	rm->pts = new point3fi[rm->wid*rm->len];
	rm->idx = new point2i[rm->wid*rm->len];
	rm->di = new BYTE[rm->wid*rm->len];
	rm->regionID = new int[rm->wid*rm->len];
	rm->segbuf = NULL;
	rm->rMap = cvCreateImage(cvSize(rm->wid, rm->len), IPL_DEPTH_8U, 3);
	rm->lMap = cvCreateImage(cvSize(rm->wid, rm->len), IPL_DEPTH_8U, 3);
}

void ReleaseRmap(RMAP *rm)
{
	delete[]rm->pts;
	delete[]rm->idx;
	delete[]rm->regionID;
	cvReleaseImage(&rm->rMap);
	cvReleaseImage(&rm->lMap);
}


