
Scene Classification  
====================

This project is an implementation of Professor James Hays's Scene Recognition with Bag of Words assignment in his (http://cs.brown.edu/courses/cs143/ "Fall 2013 class). This implementation uses multiple feature extraction techniques to classify 15 categories of sets of 100 images (total 1500 images).

Special Professor Hays for making this class's notes and assignments available open source. And also, Chun-Che Wang, Patsorn Sangkloy, Junzhe Xu, Michael Wang for posting their results on their assignments to help guide me in the right direction as I implemented my project.

Running the Project
-------------------

To run this project, you will need to take the skeleton project from the CS143 page: http://cs.brown.edu/courses/cs143/proj3/ and copy its data folder into the working directory of this project. I did not commit the images because of the amout of space it would have taken up in the git history (~87 MB). 

Some required libraries for this project are VLFeat and the Matlab's image toolbox. Once you ahve VLFeat installed you will need to replace the `run('~/Documents/MATLAB/vlfeat-0.9.19/toolbox/vl_setup')` line in proj3.m

Run `proj3` and the project will classify images in the `data/test` directory.

Spatial Pyramid Sifts
---------------------

### Finding the Pyramid of Sifts Descriptors

The first piece of creating the pyramid sift descriptors was finding the sift descriptors of the image. This was done using the dense sift implementation in VLFeat: `vl_dsift`. Creating the pyramid of sift descriptors involved splitting the image into multiple levels of sub images. That is, the first level was the full image, the second level was the image split into equal parts of four. And a third level splitting the image into 16 equal parts. Once this was done, the dense sift function was ran on the smaller portions of the images in smaller steps and appended to the image's feature vector.


### Creating the Vocabulary

For the new images to be predicted, we need to generate a vocabulary of using the sift descriptors above to create a histogram that defines the image's pyramid of sifts features.

To do this, we choose a random set of images in our data set and constructed a vocabulary using the pyramid of sifts descriptor described above. Next these descriptors were grouped together into 200 clusters using K-Means clustering with VLFeat's `vl_kmeans`. This is the most lengthy procedure in this project because of the large amounts of descriptors it needs to cluster so we have to save this vocabulary for future use in a .mat file.

### Creating the Histogram

Once the vocabulary was defined, we can define an image's pyramid of sifts with a histogram using VLFeat's `vl_kdtreequery` function. When we define an image's histogram we first run the pyramid of sift descriptors function defined initially on the image. This gives us a large vector of features defining an images sifts at multiple levels steps. We then run `vl_kdtreequery` to compare these features to our vocabularies features to give us the distance between the features we found and its closest vocabulary feature. Finally we make a histogram of these closest vocabulary features. This histogram is what we will use when we compare our training and test images.
 
Fisher Encoding of Sifts
------------------------

Fisher encoding of Sifts was similar to the pyramid of sifts descriptor, in that its goal was to preserve more spatial data of the image when running sift.

### Vocabulary

Like Pyramid Sifts, a vocabulary was generated on a random sample set of our training data. This time using mean, covariances and priors with VLFeat's gaussian mixture model implementation: `vl_gmm`.

### Encoding

Once the vocabulary was generated, we can call `vl_dsift` on our image, run it through `vl_fisher` with our means, covariances, and priors to get our fisher encoding of sifts of our image.

Gist
----

The Gist descriptor is a descriptor that describes an image by its overall spatial structure and naturalness. The implementation we use was created by Aude Oliva and Antonio Torralba here: http://people.csail.mit.edu/torralba/code/spatialenvelope/

The usage of `LMGist` was straightforward. The images were resized to a given size and a number of blocks were set to compute the gist over. Then we normalize our gist descriptors and this becomes the image's gist features.


Support Vector Machines
-----------------------

Once we have features defined for our images, we need a way to classify our test images compared to our training images. The way we use here is with support vector machines. 

### One vs Many


For each category we are trying to define, we perform a one vs many strategy against our test images. We call `vl_svmtrain` against training images' features to generate weight and offset vectors. We then multiply our weights against our test image features and add the offsets to get our one vs many result. We place this result into a prediction matrix as we go through each category and take the max to figure out which category the image is closest to our test image.


Results
-------

The results of my latest run was 80.0% using these descriptors. I have also similar results with just Fisher encoding and Gist descriptors. My assumption is that Fisher encoding and pyramid dense sifts preserves similar spatial information about an image, so using both does not provide an advantage. Below is the confusion matrix of test image categories against training image categories.

To confirm these results, I still need to perform random splits on my training and testing data and take the average, but I am confident this accuracy will remain. At least with this set of images and categories.

![Confusion Matrix](https://raw.githubusercontent.com/trngt/Scene-Classification/master/code/results_webpage/confusion_matrix.png "Confusion Matrix")


[Additional Results](https://trngt.github.io/Scene-Classification) 
 
