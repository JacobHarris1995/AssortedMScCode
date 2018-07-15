import math
import matplotlib.pyplot as plt

temp = 15.0
targTemp = 21.0 #Adjust this
htrPwr = 0.0 #the heater power

tempRecord = []

energy = 240.0
avrExtEngy = 240.0

relHum = 0.5
targRelHum = 0.7 #Adjust this

humRecord = []

def MaxHum(T): #Opperation to simulate physics, provided by group
    return 0.0018*(T/5+6.9)**3

hum = 1.0
#targHum = 1.4
humPwr = 0.0

time = 0.0

def AdjHum(temp,relHum): 
    return MaxHum(temp)*relHum

def AdjRelHum(temp,hum):
    return hum/MaxHum(temp)

def AdjTemp(energy,hum):
    return 0.1*energy/(hum+0.6)

def EnergyLoss(energy,avrExtEngy):
    return energy+(avrExtEngy-energy)**3/5000000

def AdjustHeater(temp,targTemp,htrPwr):
    if (temp < targTemp) and htrPwr < 5.0:
        return htrPwr + 1.0
    elif (temp > targTemp) and htrPwr > 0.0:
        return htrPwr - 1.0
    else:
        return htrPwr



def YananAlgorithm(temp,targTemp,htrPwr,tempRecord,time): #Yanan's Algorithm
    u = 0
    h = 1
    #a = [-3,1,2,3,4,5]
    c = [1.5,1.0,0.1]
    currentError = targTemp- temp
    PreviousError = []
    for i in tempRecord:
        PreviousError.append(targTemp-i)
    #T = (1)*(time-h)**2 + tempRecord[0] + (1)*h**2
    u = c[0]*currentError + c[1]*sum(PreviousError)/(len(PreviousError))+c[2]*(currentError - PreviousError[-1])
    u = math.atan(u)/(math.pi/2)*5
    return float(round(u))
        


def AdjustHumidifier(relHum,targRelHum,humPwr):
    if abs(relHum - targRelHum)<0.001:
        return 0.0
    elif (relHum < targRelHum) and humPwr < 2.0:
        return humPwr + 1.0
    elif (relHum > targRelHum) and humPwr > -2.0:
        return humPwr - 1.0
    else:
        return humPwr

def YananAlgorithm2(hum,temp,targRelHum,humPwr,humRecord,time): #Yanans Improved algorithm
    u = 0
    h = 1
    #a = [-3,1,2,3,4,5]
    c = [1,0.5,0.1]
    currentError = targRelHum*MaxHum(temp) - hum
    PreviousError = []
    for i in humRecord:
        PreviousError.append(targRelHum*MaxHum(temp)-i)
    #T = (1)*(time-h)**2 + tempRecord[0] + (1)*h**2
    u = c[0]*currentError + c[1]*sum(PreviousError)/(len(PreviousError))+c[2]*(currentError - PreviousError[-1])
    print(u)
    u = math.atan(u)/(math.pi/2)*5-2.5
    return float(round(u))

relHum = AdjRelHum(temp,hum)
tempRecord.append(temp)
humRecord.append(hum)

for i in range(0,250): #time progression
    print(temp,energy,relHum,'---',htrPwr,humPwr)
    time+= 1.0
    
    #htrPwr = AdjustHeater(temp,targTemp,htrPwr) #implement algorithm (temp)
    htrPwr = YananAlgorithm(temp,targTemp,htrPwr,tempRecord,time)    
    
    energy += 5.0*htrPwr
    energy = EnergyLoss(energy,avrExtEngy) #dissipation of heat
    
    temp = AdjTemp(energy,hum)
    relHum = AdjRelHum(temp,hum)
    
    if relHum >= 1: #condensation
        hum = MaxHum(temp)
        relHum = AdjRelHum(temp,hum)
    elif relHum <= 0:
        print('Error: Relative Humidity cannot be less than zero')
        break
    
    humPwr = AdjustHumidifier(relHum,targRelHum,humPwr) #implement algorithm (hum)
    #humPwr = YananAlgorithm2(hum,temp,targRelHum,humPwr,humRecord,time)    
    
    if relHum < 1 and relHum > 0:
        hum += humPwr*0.01
    relHum = AdjRelHum(temp,hum)

    if relHum > 0.6:
        if relHum > 0.7:
            hum -= (hum-MaxHum(temp)/2)**3/500
        else:
            hum -= (hum-MaxHum(temp)/2)**3/1000
    elif relHum <0.4:
        if relHum < 0.3:
            hum += (hum-MaxHum(temp)/2)**3/500
        else:
            hum += (hum-MaxHum(temp)/2)**3/1000
    tempRecord.append(temp)
    humRecord.append(hum)
plt.plot(tempRecord)
print('ping')
print(relHum)

#Write Data to file for plotting

with open("Temp.csv",mode="wt", encoding="utf8") as file: 
    for i in tempRecord:
        file.write(str(i))
        file.write(',')
    file.write(str(tempRecord[-1]))
    file.close()
with open("Hum.csv",mode="wt", encoding="utf8") as file: 
    for i in humRecord:
        file.write(str(i))
        file.write(',')
    file.write(str(humRecord[-1]))
    file.close()
print('ping')   

        
        








    
