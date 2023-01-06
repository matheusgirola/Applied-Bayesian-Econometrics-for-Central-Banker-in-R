function datapackr=packr(DATAEVIEWS)

%remove NaNs
nanindex=isnan(DATAEVIEWS);
nanindex1=sum(nanindex,2)>0;
datapackr=DATAEVIEWS;
datapackr(nanindex1,:)=[];
