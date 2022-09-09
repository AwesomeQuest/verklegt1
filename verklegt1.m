%% hluti 1
h = [31 36 32 32 31.5 33 35 31 35 36 30].*1e-2;
avgh = sum(h)/length(h);
dh = 0.5e-2;
l = 122e-2;
dl = 0.5e-2;

mu1 = avgh/l;
dmu = (((1/l)*dh)^2+(-avgh/(l^2)*dl)^2)^0.5;

%% hluti 2


data{3,2} = [];
for i = 1:3 
	for j = 1:2
		data{i,j} = readmatrix("data" + i + "h" + j);
	end
end
data{1,1} = data{1,1}(1:(end-35),:);


hold on;
for i = 1:3
	for j = 1:2
		% plot(data{i,j}(:,3))
	end
end
hold off;

samples = [];

samples{1,1} = {83:145 155:217 234:280 290:333 352:383 392:429 451:475 479:507 528:552 557:577};
samples{2,1} = {76:135 141:206 225:274 283:337 359:395 402:448 461:498 505:541 563:586 593:618};
samples{3,1} = {90:149 155:217 240:282 290:337 359:390 398:431 454:479 485:516 538:561 567:591};
samples{1,2} = {85:128 134:184 197:232 237:268 285:309 314:341 359:381 385:408 428:446 451:470};
samples{2,2} = {61:101 107:149 168:200 208:242 262:289 293:322 343:364 370:395 414:435 440:461};
samples{3,2} = {65:105 111:145 176:200 207:237 255:277 284:304 327:346 352:370 394:408 415:428};


polynoms{3,2,10} = [];

for i = 1:3
	for j = 1:2
		for k = 1:10
			ind = cell2mat(samples{i,j}(k));
			polynoms{i,j,k} = polyfit(data{i,j}(ind,1), data{i,j}(ind,3), 1);
		end
	end
end

h2 = [4.5 6.5]*1e-2;
colors = ['r' 'g'; 'b' 'c'; 'm' 'k'];
plotsvar{3,2} = [];
subplotcount = 1;
for i = 1:3
	for j = 1:2
		subplot(2,3,subplotcount)
		hold on;
		xlabel("Timi [s]")
		ylabel("Hraði [m/s]")
		title("Tilraun " + i + " við hæð " + h2(j) + " m á " + l + " m langri braut")
		plot(data{i,j}(:,1),data{i,j}(:,3), 'Color', colors(i,j))
		for k = 1:10
			ind = cell2mat(samples{i,j}(k));
			x = data{i,j}(ind,1);
			plot(x,polyval(polynoms{i,j,k}, x),'LineWidth',2, 'Color', colors(i,j))
		end
		legend("Gögn", "Línur notaðar fyrir \mu")
		hold off;
		subplotcount = subplotcount + 1;
	end
end

slopes = ones(3,2,10);

for i = 1:3
	for j = 1:2
		for k = 1:10
			temp = polynoms{i,j,k};
			slopes(i,j,k) = temp(2);
		end
	end
end


mu = ones(3,2,5);

g = 9.8;
theta = [asin(h2(1)/l) asin(h2(2)/l)];

for i = 1:3
	for j = 1:2
		for k = 1:5
			mu(i,j,k) = (slopes(i,j,k+1)-slopes(i,j,k))/(2*g*cos(theta(j)));
		end
	end
end


avgmu = sum(mu(:))/numel(mu);
maxmu = max(mu(:));
minmu = min(mu(:));

middlemu = (maxmu+minmu)/2;
dmiddlemu = middlemu-minmu;