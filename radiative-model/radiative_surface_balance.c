//Program calculates the change in surface temperature by considering the 
//radiative balance at the earth's surface. 
//
//
//Input parameters are surface temperature, dew point temperature and average 
//temperature near and under the Boundary layer inversion.
//
//Until assimilation scripts are active only surface and dewpoint temps are 
//used as input and a 10 degree F temperature difference is assumed between
//the surface and average temp near inversion height. Furthermore the pressure
//at the surface and near the inversion layer is assumed to be 1013mb and 950mb 
//respectively 
//
//Input will be normally taken from an atmospheric sounding
//
//
//Output is the surface temperature before sunrise and the total night time
//change in temperature.
//
//
//Wind and condesation or evaopration are not taken into account, however water
//vapor in the air is accounted for in the radiative and latent heat fluxes

#include <stdio.h>
#include <math.h>

#define SIG 0.000000056703     //stefan-boltzmann constant
#define EMMSFC 1               //absorption-emission constant for surface
#define TIMESTEP .1          //amount of time(seconds) that passes before new flux is calculated

//global variables
double deep_grnd_temp;     //temperature of ground under top soil surface layer
double emm;                //abosrption-emission constant for air
double fsfc;               //upward IR flux from surface
double fly1;               //downard IR flux from atmosphere
double sfctemp;            //temperature at surface
double airtemp;            //average temperature below inversion layer near BL top                       
double pottemp;            //potential temperature
double pressure1 = 101300; //surface pressure
double pressure2 = 95000;  //pressure at near inversion layer
double mixratio;           //mixing ratio
double wind = 0;           //wind speed
double dewtemp;            //dew point temperature
double beta = .4;          //moisture availability
double fract;              //fraction of clouds that cover the sky     

//functions
double calcfluxs(double sfctemp, double airtemp);         //calculates upward IR from surface
double calcfluxa(double airtemp, double sfctemp);         //calculates downward IR from atmosphere
double deltsfctemp(double flux);                         //calculates change in temperature for given timestep
double sensibleflux(double sfctemp, double airtemp);      //calculates sensible heat flux
double potential_temperature(double temp, double pressure1, double pressure2);    //calculates potential temperature
double latentflux(double sfctemp, double dewtemp);        //calculates latent heat flux
double calc_q(double rv);                                //calculates specific humidity
double calc_gflux(double sfctemp);                       //calculates flux between top soil layer and underlying earth
double mph2mpers(double wind);                           //converts miles per hour to meters per second
double Lconst(double sfctemp);                           //calculates latent heat constant taken from Satellite estimates of wind speed and latent heat flux over the global oceans., Bentamy et al.
double emission(double dtemp);                           //sets emission/absorption constant depending on water vapor content
double fahr_to_cel(double temp);                         //converts Fahrenheit to Celsius 
double cel_to_fahr(double temp);                         //converts Celsius to Fahrenheit
void calcmixingr(double dtemp);                         //calculates mixing ratio
double cloud(double fract);                              //cloud radiative parameterization                             

int main() {
    
    FILE *ofp;             //define output file

    int i;                 //loop integer 
    int time;              //amount of loops
    double sheat;           //sensible heat flux
    double latentheat;      //latent heat flux
    double rh;              //relative humidity
    double groundflux;      //flux from conduction of deep ground layer in contact with top soil
    double intemp;          //memory variable for final temperature change calculation 
    double indtemp;         //memory variable for final dew point and surface temp comparison
    double x;               //memory varialbe for relative humidity calculation
    
    ofp = fopen("out.txt", "w");

    sfctemp = 46;
    dewtemp = 30;
    fract = 0;

    printf("\nInitial Temperature = %.2f\n",sfctemp);    //prints input surface temperature
    
    airtemp = sfctemp - 11;                            //assumes 10 degree difference between surface and average temp near upper BL
    intemp = sfctemp;                                  //remembers initial surface temp to calculate total change in surface temperature  
    indtemp = dewtemp;                                 //remembers initial dew point temp to ensure realistic temperature drop

    deep_grnd_temp = sfctemp - 10;     //set underlying ground layer temperature 

    emm = emission(dewtemp);           //accounts for radiative effects of water vapor      

    if (fract > 0)
       emm = cloud(fract);

  //CONVERTS FROM FAHRENHEIT TO CELSIUS
    sfctemp = fahr_to_cel(sfctemp);
    airtemp = fahr_to_cel(airtemp);
    dewtemp = fahr_to_cel(dewtemp);
    deep_grnd_temp = fahr_to_cel(deep_grnd_temp);
     
  //CONVERTS FROM CELSIUS TO KELVIN
    sfctemp +=273;
    airtemp +=273;
    dewtemp +=273;
    deep_grnd_temp +=273;
    
   //CALCULATES INITIAL RELATIVE HUMIDITY
    calcmixingr(dewtemp);
    x = mixratio;
    calcmixingr(sfctemp); 
    rh = (x/mixratio)*100;

    printf("Initial RH = %.1f percent\n\n",rh);   //prints initial relative humidity

    time = 36000*15;                                   //sets amount of timesteps to run model

    for (i=0;i<time;i++)
    { 
      fsfc = calcfluxs(sfctemp,airtemp);            //calculates surface net radiative flux 
      sheat = sensibleflux(sfctemp,airtemp);        //calculate sensible heat flux
      latentheat = latentflux(sfctemp,dewtemp);     //calculates latent heat flux
      groundflux = calc_gflux(sfctemp);             //calculates flux from earth below surface soil layer by conduction
      fsfc = fsfc + latentheat + sheat + groundflux;//adds radiative, sensible heat, latent heat, and ground heat flux yielding net flux
      if (i%36000 == 0) {
      	printf("%.2f\n",(9*(sfctemp-273)/5)+32);
        fprintf(ofp, "%lf\n", (9*(sfctemp-273)/5)+32 );
      }
      sfctemp = sfctemp + deltsfctemp(fsfc);        //calculates change in surface temperature

    } 
      
    sfctemp -=273;                  //converts from kelvin to celsius
    sfctemp = cel_to_fahr(sfctemp);    //converts from celsius to fahrenheit
        
 
    //PRINTS FINAL SURFACE TEMP, CHANGE IN TEMP, AND HOURS OF NIGHT
    printf("\nsurface temperature at sunrise = %.2f\n", sfctemp);
    printf("Change in surface temperature  = %.2f degrees F after %d hours\n\n", sfctemp-intemp,time/36000);
 
    fclose(ofp);
    system("pause");   
    return 0;
}
/*****************************end main program********************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double calcfluxs(double sfctemp, double airtemp)
{
    double flux;       //net radiative surface flux
    
        
    flux = SIG*((emm*pow(airtemp,4)) - (EMMSFC*pow(sfctemp,4)));   //calculates flux using Stefan-Boltzmann relation
    return flux;  
}
double calcfluxa(double sfctemp, double airtemp)   //this function is not currently called upon
{
    double flux;       //net radiative flux for atmosphere near inversion layer
    
    
    flux = SIG*((EMMSFC*emm*pow(sfctemp,4)) - emm*(pow(airtemp,4)));     //calculates flux usinge Stefan-Boltzmann relation
    return flux;
}
double deltsfctemp(double flux)
{
    double newtemp;          //temperature after time step
    double csoil = 2000000;  //heat constant for layer
    double dzlay = 0.08;      //thickness of layer
    
    newtemp = (flux/(csoil*dzlay))*TIMESTEP;     //calculates change in temperature for given flux and temp step  
    return newtemp;
} 
double sensibleflux(double sfctemp, double airtemp)
{
      double density = 1; //air density
      double Cp = 1005;   //heat capicity for dry air
      double wndmix;      //temperature change from wind mixing: wind*Ch
      double sfcthta;     //surface potential temperature
      double avgthta;     //air potential temperature
      double sheat;       //snesible heat flux
      
      wndmix = 0.0025 + 0.0042*wind;                               //regression equation valid for neutral and stable BL
      sheat = density*Cp*wndmix*(airtemp - sfctemp);               //calculates sensible heat flux        

//no longer use PT to calculate sheat// avgthta = potential_temperature(airtemp,pressure1,pressure2);//calculates potential temperature      

      return sheat;
}
double latentflux(double sfctemp, double dewtemp)
{
      double density = 1;         //density of dry air
      double q;                   //actual specific humitity
      double qs;                  //saturation specific humidity
      double wndmix;              //temperature change from wind mixing: wind*Ch
      double latentheat;          //total latent heat flux
      double lhcnst;              //latent heat of vaporization constant = 2501000 J/kg at 0c
                                 //latent heat of saturation const = 2834000 J/kg
                                 //latent heat of fusion const = 333700 J/kg

      wind = mph2mpers(wind);            //converts wind from mph to meters per second
      wndmix = 0.0025 + 0.0042*wind;    //regression equation valid for neutral BL
      lhcnst = Lconst(sfctemp);          //calculates latent heat of evaporation
      calcmixingr(sfctemp);              //calculates saturation mixing ratio
      qs = calc_q(mixratio);             //calculates saturation specific humidty
      calcmixingr(dewtemp);              //calculates mixing ratio
      q = calc_q(mixratio);              //calculates specific humidty
      
      latentheat = density*wndmix*beta*lhcnst*(q - qs); //calculates latent heat flux
        
      return latentheat;
}         
double potential_temperature(double temp, double pressure1, double pressure2)
{
     double kdry;    //poisson constant for dry atmosphere
     double kmoist;  //poisson constant for moist atmosphere
     double pottemp; //potential temperature
     double pavg;    //average atmospheric pressure
     
   //INITIALIZE POISSON CONSTANT
     kdry = 0.2854;                          
     kmoist = 0.2854*(1 - 0.24*mixratio);
     
     pavg = ((0.7*pressure1)+pressure2)/2;        //calculates simple average press     
     pottemp = temp*(pow((pressure1/pavg),kdry)); //calculates potential temperature     
     return pottemp;
}
void calcmixingr(double dtemp)
{
     double e;     //vapor pressure

     dtemp = dtemp - 273;                              //converts from Kelvin to Celsuis
     e = 6.11*(pow(10,((7.5*dtemp)/(237.7+dtemp))));   //converts from dew point temp to vapor pressure
     e = e*100;                                        //converts from hPa to Pa
     mixratio = (0.622*e)/(pressure1 - e);             //computes mixing ratio 
     mixratio = mixratio*1;                            //convert to g/Kg   
}
double calc_q(double rv)
{
     double specific_humidity;              //define specific humidity variable
     specific_humidity = rv/(1 + rv);      //calculates specific humidity
     return specific_humidity;
}
double calc_gflux(double sfctemp)
{
      double k = 1;                //thermal conductivity parameter
      double dz;               //depth of layer between soil surface and deep soil layer
      double Gflux;            //flux in watts per meter^2 between top and bottom soil layers

          
      Gflux = (k*(deep_grnd_temp - sfctemp)/0.8);   //calculates flux from deep ground layer
      return Gflux;
}
double emission(double dtemp)
{
    double emm;
    
    //ASIGNS RADIATIVE PARAMETER BASED ON MOISTURE IN THE AIR
    //MEASURED BY THE DEW POINT TEMPERATURE
    if (dewtemp >= 85)
       emm = 0.92;             
    else if (85 > dewtemp && dewtemp >= 70)
         emm = 0.88;         
    else if (70 > dewtemp && dewtemp >= 50)
         emm = 0.85;    
    else if (50 > dewtemp && dewtemp >= 35)
         emm = 0.75;         
    else if (35 > dewtemp && dewtemp >= 15)
         emm = 0.71;
    else if (15 > dewtemp && dewtemp >= 0)
         emm = 0.64;
    else if (0 > dewtemp && dewtemp >= -15)
         emm = 0.49;
    else if (-15 > dewtemp )
         emm = 0.45;
    return emm;         
}
double cloud(double fract)
{
    //MODIFIES RADATIVE BALANCE DEPENDING ON CLOUD COVER
    if (fract >= 0.9)
       emm = 1;             
    else if (0.9 > fract && fract >= 0.8)
         emm = 0.9;         
    else if (0.8 > fract && fract >= 0.7)
         emm = 0.85;         
    else if (0.7 > fract && fract >= 0.6)
         emm = 0.75;         
    else if (0.6 > fract && fract >= 0.5)
         emm = 0.65;         
    else if (0.4 > fract && fract >= 0.3)
         emm = emm*1.086956;   
    return emm;
}  
double Lconst(double sfctemp)
{
      double Lheat;
      sfctemp -=273;                              //converts from kelvin to celsius
      Lheat = 4186.8*(597.31 - 0.5625*sfctemp);   //calculates latent heat constant
      return Lheat;
}
double mph2mpers(double wind)
{
     wind = ((wind*1.6*1000)/3600);                 //converts wind from mph to meters per second
     return wind;
}  
double fahr_to_cel(double temp)
{   
   temp = (5*(temp-32))/9;
   return temp;
}
double cel_to_fahr(double temp)
{
   temp = ((temp*9)/5) + 32; //converts from celsuis to farhrenheit
   return temp;
}

