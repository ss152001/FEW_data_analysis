def get_current_discharge():
    # retrieve current discharge data of all stations from USGS
    time = []
    current_Q = []
    import urllib
    import pandas as pd 
    import numpy as np
    import sys,os
    #setup path
    dirname = "FEW_data_analysis"
    datapath = os.getcwd()
    path = os.path.join(datapath, dirname)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    filename = '/stations.csv'
    d = pd.read_csv(path+filename)
    for i in range(len(d)):
        site_number = d['site_number'][i]
        if len(str(site_number)) ==7:
            site_number = '0'+str(site_number)
        state_code = d['code'][i]
        historical_days = 1
        url='https://nwis.waterdata.usgs.gov/'+state_code+'/nwis/uv?cb_00060=on&format=rdb&site_no='+str(site_number)+'&period='+str(historical_days)
        filename = str(site_number)+'.txt'
        urllib.request.urlretrieve(url,path+filename)
        f = urllib.request.urlretrieve(url)
        df = pd.read_csv(path+filename,sep='\t',comment='#')
        for j in range(len(df.columns)):
            c = df.columns[j]
            if c[-5:] == '00060':
                index = j
        Q = [_ for _ in df[df.columns[index]]]
        Q = Q[-1]
        T = [_ for _ in df['datetime']]
        T = T[-1]
        time.append(T)
        number = [str(_) for _ in np.arange(10)]
        if Q[0] in number:
            current_Q.append(round(float(Q)*0.0283168,2))
        else:
            Q == 0
            current_Q.append(round(float(Q)*0.0283168,2))
        print(round(i*100/len(d),2),'% completed')
    time = pd.DataFrame(time)
    time.columns = ['datetime']
    current_Q = pd.DataFrame(current_Q)
    current_Q.columns = ['Q(m3s-1)']
    data = pd.concat([d,time,current_Q],axis=1)
    data.to_csv(path+'/current_discharge.csv')
    return data


    





      
        

        
def get_elevation(long,lat):
    # input: long, lat
    # return: elvation (m)
    import requests
    import urllib
    import pandas as pd
    url = r'https://nationalmap.gov/epqs/pqs.php?'
    elevation=[]
    params = {
        'output': 'json',
        'x': long,
        'y': lat,
        'units': 'Meters'
    }
    result = requests.get((url + urllib.parse.urlencode(params)))
    elevation.append(result.json()['USGS_Elevation_Point_Query_Service']['Elevation_Query']['Elevation'])
    z = elevation[0]
    return z

def get_distance(long1,lat1,long2,lat2):
    # input: long, lat of 2 locations
    # return distance (m)
    from geopy import distance
    coords_1 = (lat1,long1)
    coords_2 = (lat2,long2)
    distance = distance.distance(coords_1, coords_2).m
    return distance



def get_L_dH(site_number):
    # input: site number (e.g., 12340500)
    # L: stream length
    # return: stream length (m) and gross head (m)
    # plot 3D local topography
    import os
    dirname = "FEW_data_analysis"
    datapath = os.getcwd()
    path = os.path.join(datapath, dirname)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

    import pandas as pd 
    import numpy as np
    filename = '/stations.csv' 
    d = pd.read_csv(path+filename)
    for i in range(len(d)):
        if str(site_number) == str(d['site_number'][i]):
            lat = d['lat'][i]
            long = d['long'][i]
            h,k = long,lat
            r = 0.0025 # set stream length at 500m
            X = []
            Y = []  
            Z = []          
            angle = np.arange(0,360,1)
            for a in angle:
                x = r*np.cos(np.radians(a))
                y = r*np.sin(np.radians(a))
                X.append(x+h)
                Y.append(y+k)
                z = get_elevation(x+h,y+k)
                Z.append(z)
                print('calculating gross head...',round(a*100/360,2),'%')
            z = pd.concat([pd.DataFrame(X),pd.DataFrame(Y),pd.DataFrame(Z)],axis=1)
            z.columns = ['long','lat','Z']      
            for i in range(len(z)):
                if z['Z'][i] == min([j for j in z['Z']]):
                    x_down,y_down,z_down = z.iloc[i][0],z.iloc[i][1],z.iloc[i][2]
            def get_m(x1,y1,x2,y2):
                dy = y2-y1
                dx = x2-x1
                if x2-x1 == 0:
                    dx = dx+ 0.001
                m = dy/dx
                return m                      
            m1 =get_m(h,k,x_down,y_down)
            collect_m2 =[]
            m1m2 = []
            i_ = []
            for i in range(len(z)):
                xx1 = z['long'][i]
                yy1 = z['lat'][i]
                m2 = get_m(h,k,xx1,yy1)
                collect_m2.append(m2)
                m1m2_ = round(m1*m2,2)
                m1m2.append(m1m2)
                if m1m2_ == -1.0:
                    i_.append(i)
            z_upstream = [i for i in z['Z'][i_[0]:i_[1]]] 
            x_upstream = [i for i in z['long'][i_[0]:i_[1]]] 
            y_upstream = [i for i in z['lat'][i_[0]:i_[1]]] 
            
            up = pd.concat([pd.DataFrame(x_upstream),pd.DataFrame(y_upstream),pd.DataFrame(z_upstream)],axis=1)
            up.columns = ['long','lat','Z']    
            for i in range(len(up)):
                min_z = min(up['Z'])
                if up['Z'][i] == min_z:
                    x_up,y_up,z_up = up['long'][i],up['lat'][i],up['Z'][i]  
            
            from mpl_toolkits import mplot3d
            import matplotlib.pyplot as plt
            import matplotlib as mpl
            from mpl_toolkits.mplot3d import Axes3D
            z_ = np.array(Z)
            x_ = np.array(X)
            y_ = np.array(Y)
            ax = plt.axes(projection='3d')
            zdata = z_
            xdata = x_
            ydata = y_
            L = get_distance(x_down,y_down,x_up,y_up)
            dH = get_elevation(x_up,y_up)-get_elevation(x_down,y_down)
            L = round(L,4) #stream length (m)
            dH = round(dH,6) #gross headh (m)
            if L == 0:
                z_upstream = [i for i in z['Z'][0:i_[0]]] + [i for i in z['Z'][i_[1]+1:]]
                x_upstream = [i for i in z['long'][0:i_[0]]] + [i for i in z['long'][i_[1]+1:]]
                y_upstream = [i for i in z['lat'][0:i_[0]]] + [i for i in z['lat'][i_[1]+1:]]
                up = pd.concat([pd.DataFrame(x_upstream),pd.DataFrame(y_upstream),pd.DataFrame(z_upstream)],axis=1)
                up.columns = ['long','lat','Z']    
                for i in range(len(up)):
                    min_z = min(up['Z'])
                    if up['Z'][i] == min_z:
                        x_up,y_up,z_up = up['long'][i],up['lat'][i],up['Z'][i]  
                z_ = np.array(Z)
                x_ = np.array(X)
                y_ = np.array(Y)
                ax = plt.axes(projection='3d')
                zdata = z_
                xdata = x_
                ydata = y_
            ax.scatter3D(xdata, ydata, zdata, s=5,c=zdata, cmap='Greens')
            ax.scatter3D(x_down,y_down,z_down,c='r')
            ax.text(x_down,y_down,z_down,'downstream',fontsize=8)
            ax.scatter3D(h,k,get_elevation(h,k),c='r',marker ='*')
            ax.text(h,k,get_elevation(h,k),'gaging station',fontsize=8)
            ax.scatter3D(x_up,y_up,z_up,c='r')
            ax.text(x_up,y_up,z_up,'upstream',fontsize=8)
            ax.plot3D([x_down,x_up], [y_down,y_up], [z_down,z_up], 'b')
            ax.set_title('Local topography at station: '+str(site_number))
            ax.set_xlabel('long')
            ax.set_ylabel('lat')
            ax.set_zlabel('elevation(m)')
            ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.3f}'))
            ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.3f}'))
            ax.zaxis.set_tick_params(labelsize=10)
            ax.xaxis.set_tick_params(labelsize=5)
            ax.yaxis.set_tick_params(labelsize=5)
            ax.view_init(20,azim=60)
            L = get_distance(x_down,y_down,x_up,y_up)
            dH = get_elevation(x_up,y_up)-get_elevation(x_down,y_down)
            L = round(L,4) #stream length (m)
            dH = round(dH,6) #gross headh (m)
            print('stream length: ',L,' m','gross head: ',dH,' m')
    return L,dH


def plot_historical_discharge_energy(site_number,year1,year2):
    # input: site number
    # year1: starting year
    # year2: ending year
    stream_length,dH = get_L_dH(site_number) #get gross head (m)
    # retrieve current discharge data of all stations from USGS
    import urllib
    import pandas as pd 
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    #setup path
    dirname = "FEW_data_analysis"
    datapath = os.getcwd()
    path = os.path.join(datapath, dirname)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
    filename = '/stations.csv'
    d = pd.read_csv(path+filename)
    if len(str(site_number)) ==7:
        site_number = '0'+str(site_number)
    for i in range(len(d)):
        if d['site_number'][i] == site_number:
            state_code = d['code'][i]
            historical_days = ((year2-year1)+3)*365
    url='https://nwis.waterdata.usgs.gov/'+state_code+'/nwis/uv?cb_00060=on&format=rdb&site_no='+str(site_number)+'&period='+str(historical_days)
    filename = str(site_number)+'.txt'
    urllib.request.urlretrieve(url,path+filename)
    f = urllib.request.urlretrieve(url)
    df = pd.read_csv(path+filename,sep='\t',comment='#')
    for j in range(len(df.columns)):
            c = df.columns[j]
            if c[-5:] == '00060':
                index = j
    T = [_ for _ in df['datetime']]
    T = T[1:]
    Q = [_ for _ in df[df.columns[index]]]
    Q = Q[1:]
    q = []
    number = [str(_) for _ in np.arange(10)]
    for i in range(len(Q)):
        Qi=str(Q[i])
        if Qi[0] in number:
            q.append(round(float(Q[i])*0.0283168,2)) #convert cfs to m3s-1
        else:
            q.append(0)   
    Tq = pd.concat([pd.DataFrame(T),pd.DataFrame(q)],axis=1)
    Tq.columns = ['time','Q(m3s-1)']
    d_ = {}
    for j in range(len(Tq)):
        k = Tq['time'][j][0:10] 
        v = Tq['Q(m3s-1)'][j]
        d_[k] = v
    fig = plt.figure()
    ax = fig.add_subplot(111)
    Qmax = []
    Emax = []
    for i in np.arange(year1,year2+1,1):
        year_i = str(i)
        Y = []
        E =[]
        #convert gross head to energy
        for k,v in d_.items():
            if k[0:4] == year_i:
                Y.append(v)
        E = [9.8*dH*0.9*Q for Q in Y]
        import matplotlib.pyplot as plt
        Emax.append(max(E))
        Qmax.append(max(Y))
        ax.plot(Y,label=year_i)
        ax.set_xlim(0,365)
        ax.set_ylim(0,max(Qmax)+20)
        ax.set_xlabel('Day of year')
        ax.set_ylabel(r'Q ($m^{3} s^{-1}$)')
        ax.legend(fancybox=True,framealpha=0)
    ax2 = ax.twinx()
    ax2.grid()
    ax2.set_ylabel('Power potential (kW)')
    plt.title('Data from '+str(year1)+' - '+ str(year2)+' of site number:'+str(site_number))
    ax2.set_ylim(0,max(Emax))      

    
plot_historical_discharge_energy(12340500,2010,2018)

plot_historical_discharge_energy(12388700,2012,2018)
plot_historical_discharge_energy(7022000,2012,2018)
