T10N = fix(10*fs/f1estim);
at = 0.1*fs+T10N;


[precval,precidx] = max(-comp0(1:at));

figure(1); clf; plot(pipecurr); xlim([1 1e4]);grid on;
hold on; plot(at*[1,1],20*[-1,1],'--r');
plot(precidx*[1,1],20*[-1,1],'--g');
plot( (precidx+fix(fs/f1estim))*[1,1],20*[-1,1],'--g');
plot( (precidx-fix(fs/f1estim))*[1,1],20*[-1,1],'--g');



figure(2); clf; plot(comp0);grid on; xlim([1 1e4]);
hold on; plot(at*[1,1],20*[-1,1],'--r');
hold on; plot(precidx*[1,1],20*[-1,1],'--g');
plot( (precidx+fix(fs/f1estim))*[1,1],20*[-1,1],'--g');
plot( (precidx-fix(fs/f1estim))*[1,1],20*[-1,1],'--g');

figure(3); clf; plot(COMP(1,:)); grid on; xlim([1 1e4]);
hold on; plot(at*[1,1],20*[-1,1],'--r');
hold on; plot(precidx*[1,1],20*[-1,1],'--g');
plot( (precidx+fix(fs/f1estim))*[1,1],20*[-1,1],'--g');
plot( (precidx-fix(fs/f1estim))*[1,1],20*[-1,1],'--g');

figure(4); clf; plot(foot_trans{idx});xlim([1 1e4]);
hold on; plot(precidx*[1,1],20*[-1,1],'--g');
plot( (precidx+fix(fs/f1estim))*[1,1],20*[-1,1],'--g');
plot( (precidx-fix(fs/f1estim))*[1,1],20*[-1,1],'--g');
