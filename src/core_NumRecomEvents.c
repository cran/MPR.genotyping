void core_NumRecomEvents(int *h,int *Col,int *Row,int *R)
{
    int i=0,j=0,worknow=0;
    R[0]=0;
    for (i=0;i<Row[0];i++)//i����R���������
    {
        worknow=h[i*Col[0]+0];
        if (worknow==0)
        {
            R[0]-=1;//��Ϊһ�в�����һ��SNP��û�в⵽�������Կ�����˴���
        }
        for (j=1;j<Col[0];j++)//����
        {
            if(h[i*Col[0]+j]!=0)//�жϲ���0
            {
                if(h[i*Col[0]+j]!=worknow)//���Һͱ�׼ֵ����ͬ
                {
                    R[0]+=1;
                    worknow=h[i*Col[0]+j];
                }
            }
        }
        if (worknow==0)
        {
            R[0]+=1;
        }
    }
}
