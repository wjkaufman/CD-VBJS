function [] = savefig(name)

set(gcf, 'PaperPosition', [0 0 7 5]);
set(gcf, 'PaperSize', [7 5]);
print(name, '-dpdf');

end