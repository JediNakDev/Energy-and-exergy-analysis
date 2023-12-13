from powerdata import NoPCM,P1,P2,P3,P4,P5,P6,P7,M12,M13,M24,M26,M35,M36,M37
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
x = []
for i in range(150,570):
    time=datetime(2022,3,1,i//30,(i%30)*2)
    x.append(time)
dP1=[]
dP2=[]
dP3=[]
dP4=[]
dP5=[]
dP6=[]
dP7=[]
dP12=[]
for i in range(420):
    dP1.append(P1[i]-NoPCM[i])
    dP2.append(P2[i]-NoPCM[i])
    dP3.append(P3[i]-NoPCM[i])
    dP4.append(P4[i]-NoPCM[i])
    dP5.append(P5[i]-NoPCM[i])
    dP6.append(P6[i]-NoPCM[i])
    dP7.append(P7[i]-NoPCM[i])
    dP12.append(M12[i]-NoPCM[i])
ax = plt.subplot()
ax.plot_date(x, P1,'k',xdate=True)
ax.plot_date(x, NoPCM,'k-.',xdate=True)
# ax.plot_date(x, dP1,'k',xdate=True)
# ax.plot_date(x, dP2,'k',dashes=[1,1.5],xdate=True)
# ax.plot_date(x, dP3,'k',dashes=[4, 1],xdate=True)
# ax.plot_date(x, dP4,'k-.',xdate=True)
# ax.plot_date(x, dP5,'gray',xdate=True)
# ax.plot_date(x, dP6,'gray',dashes=[1,1.5],xdate=True)
# ax.plot_date(x, dP7,'gray',dashes=[4, 1],xdate=True)
ax.set_ylabel('Power area density [W/$m^{2}$]',fontsize='36')
ax.set_xlabel('Time',fontsize='36')
ax.tick_params(axis='both', which='major', labelsize=32)
myFmt = mdates.DateFormatter('%H:%M')
ax.xaxis.set_major_formatter(myFmt)
# plt.legend(['Paraffin Wax','Sodium Sulfate Decahydrate','Polyethylene Glycol','Calcium Chloride Hexahydrate','Octadecane','Lauric Acid','Capric Acid'],loc='upper right',fontsize=14)
plt.ylim(0,250)
# plt.legend(['PV-PCM','Conventional PV'],loc='upper right',fontsize='36')
plt.grid(True,alpha=0.5)
plt.show()