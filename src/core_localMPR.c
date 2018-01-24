int NumRecomEvents(int *h, int *Row, int *Col);
int *core_localMPR(int *SNPmatrix,int *Row, int *Col, int *maxNStep, int *allelelog){

    int NStep,i,j,k,s,GlobalBestcountR,NowcountR,BestcountR;
    NStep=1;
    GlobalBestcountR=Row[0]*Col[0];
    BestcountR=Row[0]*Col[0];
    while(NStep<=maxNStep[0])//�ı����
    {
        GlobalBestcountR=BestcountR;

        BestcountR=NumRecomEvents(SNPmatrix,Row,Col);
        for (i=0;i<Row[0]-NStep+1;i++)
        {
            //
            k=i;
            for (s=0;s<NStep;s++)
            {
                //����ȡ�෴��
                for (j=0;j<Col[0];j++)
                {
                    SNPmatrix[k*Col[0]+j]=-SNPmatrix[k*Col[0]+j];
                }
                k += 1;
            }
            //
            NowcountR=NumRecomEvents(SNPmatrix,Row,Col);
            //printf("%d!",NowcountR);
            if (NowcountR>=BestcountR)//�Ⲩ���þ͵ñ����
            {
                //
                k=i;
                for (s=0;s<NStep;s++)
                {
                    //����ȡ�෴��
                    for (j=0;j<Col[0];j++)
                    {
                        SNPmatrix[k*Col[0]+j]=-SNPmatrix[k*Col[0]+j];
                    }
                k += 1;
                }
            //
            }
            else//���ı�Ͳ���¼��
            {
                BestcountR=NowcountR;
                for (s=0;s<NStep;s++)
                {
                    allelelog[i+s]=-allelelog[i+s]+1;
                }
            }
        }

        //���NStep=1�ã��Ǿ�������������һ��NStep=1
        if (GlobalBestcountR>BestcountR)
        {
            NStep=1;
        }
        else
        {
            NStep++;
        }
    }

    return allelelog;
}
int NumRecomEvents(int *h, int *Row, int *Col)
{
    /*
    printf("%d\n",h[0]);
    printf("%d\n",h[1]);
    printf("%d\n",h[2]);
    */
    int i,j,Rcount=0;
    char worknow;

    for (j=0;j<Col[0];j++)
    {
        worknow=h[0*Col[0]+j];
        if (worknow==0)
        {
            Rcount-=1;//��Ϊһ�в�����һ��SNP��û�в⵽�������Կ�����˴���
        }
        //printf("%d\n",h[j]);
        for (i=1;i<Row[0];i++)//����
        {
            if(h[i*Col[0]+j]!=0)//�жϲ���0
            {
                if(h[i*Col[0]+j]!=worknow)//���Һͱ�׼ֵ����ͬ
                {
                    Rcount+=1;
                    worknow=h[i*Col[0]+j];
                    //printf("[%d-%d]",i,j);
                }
            }
        }

        if (worknow==0)
        {
            Rcount+=1;
        }
        //
    }
    //printf("%d\n",Rcount);
    return Rcount;
}
