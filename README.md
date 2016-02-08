# particle_tracking
Code to track monodisperse disks from images. 
Main tracking code is in the function complete_track(folder number, experiment number)
You give the folder and experiment where the images are stored. 
To be more efficient uses matlab 2011 parallel computer resources. 

The code divides in 4 main parts:

1) Get background. 
2) Find the centers.
3) Organize the files so the have continuity.
4) Connect the centers.
5) Identify the particles that move.


1) To get the background complete_track rans the function getbackground inside the scheduler. 
In it first run, it opens each .mat file that contains a series of image frames. For each image runs a filter using the discs size and keeps only the pixels above a cutoff value. The pixels that contains disks parts are darker, we only want the brighter pixels that are part of the background. We normalize each image to cancel changes over time in illumination. It then keeps the brightest pixels of each image, as they are the ones that are not affected by the presences of the particles. 
At the end of the first run complete_image loads the image created from each experimental file and then creates a new image bk1 from the brightest pixels of each. 

In it second run, it uses bk1 and get the image difference between the image of each file and bk1 to create bk2 that will be used to normalize the image-background for each image. 

2) To find the centers of each particle complete_track runs the function findparticles center and passes the information about the files where the images are stored as well as the diameter and intensity width of each particle and the brightness cutoff to find the peaks.  Find particles center first normalize each image using bk1 and bk2 so the new value of the pixels in the new image sim that is between 0 and 1. Where ideally the background pixels have a zero value and the disks pixels one. findparticle centers.m then runs chiimg that returns a matrix ichi of the image size with the difference of the image sim with the image of an ideal particle created by function ipf. Low values of ichi indicates that at that position the image difference is small. So the particle centers are located at the minimal points of ichi. The function findpeaks(mk./ichi,mk,Cutoff,MinSep) finds the values of the image mk./ichi that are above a cutoff. Then it saves only the positions that are at least Minsep away from each other. If two points are closer than MinSep, the function returns the position with the higher value.  mk is just a matrix of ones and zeros to get rid of parts of the image were we are not looking for particles. 
The returned positions are the center of the particles. After finding the particles the function findpeaks compares the first and last image of each file and check to see if the particles moved. 

3) complete_track uses the information from findtracks to organize the files. It sticks together consecutive files where the particles moved. 

4) Uses function mytrack. mytrack receives the information of from which file to which file connect, then it opens the files and loads the position information of each particles at each time step and find the minimal total displacement of the particles between two consecutive images. Mytrack creates the files Tracked with the position of the particles PX, PY at each time step. So PX(1,t) gives the trajectory of particle one. Mytrack also writes warning files for when the code finds a particle in frame t but not in frame t+1 and save the ghost particle position as well as the time step when this happens. If it finds no movement at all the it creates a file called no avalanches. At the end of this step the particles are tracked and information can be retrieved. 


5) Finally it calls the function displacement. Displacement loads the data in the files previously created and compares the position of the particles in the first and last frame to find out which particles moved. It then finds the velocities at each time step comparing the positions every n numbers of frames and applying a filter to the signal to smear the noise that comes from the errors in findparticles code. If the file tracked does not exists gives a warning indicating that something went wrong in the mytrack part of the code. Finally the function avalanche_size collects the data from all the files and counts the number of avalanches, size of each one etc. The results are saved in new files. One per experimental run. 





