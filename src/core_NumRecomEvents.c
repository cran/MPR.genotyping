void core_NumRecomEvents(int *h,int *Col,int *Row,int *R)
{
    int i=0,j=0,worknow=0;
    R[0]=0;
    for (i=0;i<Row[0];i++)//i代表R里面的行数
    {
        worknow=h[i*Col[0]+0];
        if (worknow==0)
        {
            R[0]-=1;//因为一列不可能一个SNP都没有测到过，所以可以如此处理
        }
        for (j=1;j<Col[0];j++)//行数
        {
            if(h[i*Col[0]+j]!=0)//判断不是0
            {
                if(h[i*Col[0]+j]!=worknow)//并且和标准值不相同
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
