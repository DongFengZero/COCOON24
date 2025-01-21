# An Approach to Tight I/O Lower Bounds for Algorithms with Composite Procedures

Due to space limitations, we provide the paper's appendix in this repository.

Our exploratory experiments on NAS tasks are particularly grateful for the work on Efficient Neural Architecture Search (ENAS).
https://github.com/carpedm20/ENAS-pytorch

For readers who are interested in "Introducing the I/O Lower Bound Theorem into NAS Tasks", we strongly recommend that you read the appendix first, where we explain the experimental setup and the suggestions and future directions we give based on our experimental experience. If you have more suggestions or questions, you are welcome to discuss them via email.

If you have further questions, please contact xiarui21@nudt.edu.cn

Taking into account the randomness of the experimental results, we provide all the neural network models generated in 100 (101) epochs of the experiment that generated the data given in the paper in TEST.zip for verification.

For more information, see 
Address：https://pan.baidu.com/s/1Z0M8JnUaqcWcoA1FCvXpYA?pwd=6nc7 
Code：6nc7 

Note: We provide the training model weights that produced the corresponding results. Similar results may be produced based on these weights (not guaranteed to be the same. For completely identical results, please refer to the 100 (CIFAR100：101) rounds of neural network structure diagram we provided). These files serve as evidence that we implemented the experiment and are for researchers' reference only (CIFAR100 is 99-100 rounds, CIFAR100 is 108-109 rounds. The weights of CIFAR100 training at 100 epochs were automatically cleared and could not be saved, so the model weights of 108-109 epochs when the training was stopped late are provided):
