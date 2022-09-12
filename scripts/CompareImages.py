from skimage.metrics import structural_similarity as ssim
import matplotlib.pyplot as plt
import numpy as np
import cv2
import sys
def mse(imageA, imageB):
	# the 'Mean Squared Error' between the two images is the
	# sum of the squared difference between the two images;
	# NOTE: the two images must have the same dimension
	err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
	err /= float(imageA.shape[0] * imageA.shape[1])

	# return the MSE, the lower the error, the more "similar"
	# the two images are
	return err
def compare_images(imageA, imageB):
    # compute the mean squared error and structural similarity
    # index for the images
    m = mse(imageA, imageB)
    #s = ssim(imageA, imageB)
    print(m)
    #setup the figure
    fig = plt.figure("Compare")
    plt.suptitle("MSE: %.2f" % (m))#, SSIM: %.2f" % (m, s))
    # show first image
    ax = fig.add_subplot(1, 3, 1)
    plt.imshow(imageA, vmin=0,vmax=255,cmap = plt.cm.gray)
    plt.axis("off")
    # show the second image
    ax = fig.add_subplot(1, 3, 2)
    plt.imshow(imageB,  vmin=0,vmax=255,cmap = plt.cm.gray)
    plt.axis("off")
    ax = fig.add_subplot(1, 3, 3)

    help = abs(imageA-imageB);

    plt.imshow(help, cmap = plt.cm.plasma)
    plt.axis("off")

    plt.show()


if(len(sys.argv) < 3):
    print("Please use the filenames as args")

im1 = cv2.imread(sys.argv[1])
im1 = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)


im2 = cv2.imread(sys.argv[2])
im2 = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)

compare_images(im1,im2)