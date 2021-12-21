clear classes;
clc;

convolutional_encoder = ConvEncoder.E2;

arrayfun(@(ebn0) sum(arrayfun(@(d) convolutional_encoder.enumerate_weights(d) * qfunc(sqrt(2*d*rate*ebn0)), scale)), 10.^(EbN0 / 10));
figure()
plot(thr())